!-----------------------------------------------------------------------
      MODULE control_flags
!-----------------------------------------------------------------------
        IMPLICIT NONE
        SAVE
        CHARACTER(LEN=80) :: job_type
        INTEGER :: nfil
        INTEGER :: lpr
        INTEGER :: immb, itp
        LOGICAL :: implicit_fluxes, implicit_enthalpy
!-----------------------------------------------------------------------
      END MODULE
!-----------------------------------------------------------------------
