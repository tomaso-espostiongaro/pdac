!-----------------------------------------------------------------------
      MODULE control_flags
!-----------------------------------------------------------------------
        IMPLICIT NONE
        SAVE
        CHARACTER(LEN=80) :: job_type
        INTEGER :: nfil
        INTEGER :: lpr
        INTEGER :: imr
        LOGICAL :: implicit_fluxes
        LOGICAL :: implicit_enthalpy
        LOGICAL :: formatted_output
        CHARACTER(LEN=4) :: prog   ! "PDAC" "POSP"
!-----------------------------------------------------------------------
      END MODULE
!-----------------------------------------------------------------------
