!-----------------------------------------------------------------------
      MODULE control_flags
!-----------------------------------------------------------------------
        IMPLICIT NONE
        SAVE
        CHARACTER(LEN=80) :: job_type
        INTEGER :: lpr
        INTEGER :: imr
        LOGICAL :: implicit_fluxes
        LOGICAL :: implicit_enthalpy
        LOGICAL :: run
!-----------------------------------------------------------------------
      END MODULE
!-----------------------------------------------------------------------
