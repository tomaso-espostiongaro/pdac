!----------------------------------------------
      MODULE roughness_module
!----------------------------------------------
      IMPLICIT NONE
      REAL*8 :: zrough
      REAL*8, ALLOCATABLE :: roughness(:)
      REAL*8, ALLOCATABLE :: roughness2d(:,:)
      SAVE
!----------------------------------------------
      CONTAINS
!----------------------------------------------
      SUBROUTINE roughness_setup
      USE control_flags, ONLY: job_type
      USE dimensions 
      IMPLICIT NONE
!
      IF (job_type == '2D') THEN
              ALLOCATE(roughness(nx))
              roughness = zrough
      ELSE IF (job_type == '3D') THEN
              ALLOCATE(roughness2d(nx,ny))
              roughness2d = zrough
      ELSE
              CALL error('roughness','unknown job type',1)
      END IF
!
      RETURN
      END SUBROUTINE roughness_setup
!----------------------------------------------
      END MODULE
!----------------------------------------------
