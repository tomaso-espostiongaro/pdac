!-----------------------------------------------------------------------
      MODULE control_flags
!-----------------------------------------------------------------------
        IMPLICIT NONE
        SAVE
        CHARACTER(LEN=80) :: job_type
        INTEGER :: nfil = 0
        INTEGER :: lpr
        LOGICAL :: implicit_fluxes, implicit_enthalpy
!-----------------------------------------------------------------------
      END MODULE
!-----------------------------------------------------------------------
