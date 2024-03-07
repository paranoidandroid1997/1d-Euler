module sim_data

#include "definition.h"

    implicit none

  !! numerics
    real, save :: sim_cfl, sim_tmax, sim_outputIntervalTime
    integer, save :: sim_order, sim_nStep
    character(len=MAX_STRING_LENGTH), save :: sim_name, sim_limiter, sim_riemann, IC_type
    logical, save :: sim_charLimiting

  !! ICs
    real, save :: sim_gamma
    real, save :: sim_densL, sim_velxL, sim_presL
    real, save :: sim_densR, sim_velxR, sim_presR
    real, save :: sim_densM, sim_velxM, sim_presM
    real, save :: sim_shockLoc, sim_shockLoc2
    real, save :: sim_smallPres

  !! BCs
    character(len=MAX_STRING_LENGTH), save :: sim_bcType

  !! IO
    integer, save :: sim_ioNfreq
    real, save :: sim_ioTfreq

end module sim_data
