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
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions 
      USE volcano_topography, ONLY: itp, topo_c
      IMPLICIT NONE
      INTEGER :: i
!
      IF (job_type == JOB_TYPE_2D) THEN
              ALLOCATE(roughness(nx))
              DO i=1,nx
              roughness(i) = zrough
!                IF (itp > 1) THEN
!                      IF (topo_c(i) == 0.D0) roughness(i) = 0.05D0
!                END IF
              END DO
      ELSE IF (job_type == JOB_TYPE_3D) THEN
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
