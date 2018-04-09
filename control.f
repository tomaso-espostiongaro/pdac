!-----------------------------------------------------------------------
      MODULE control_flags
!-----------------------------------------------------------------------
        IMPLICIT NONE
        SAVE
        INTEGER :: job_type
        INTEGER, PARAMETER :: JOB_TYPE_3D = 3
        INTEGER, PARAMETER :: JOB_TYPE_2D = 2
        INTEGER, PARAMETER :: JOB_TYPE_1D = 1
        INTEGER :: nfil
        LOGICAL :: cstop
        INTEGER :: lpr
        INTEGER :: imr
        LOGICAL :: implicit_fluxes
        LOGICAL :: implicit_enthalpy
        LOGICAL :: formatted_output
        LOGICAL :: formatted_input
        CHARACTER(LEN=4) :: prog   ! "PDAC" "POSP"
!-----------------------------------------------------------------------
      END MODULE
!-----------------------------------------------------------------------
