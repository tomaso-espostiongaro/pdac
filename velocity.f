!-----------------------------------------------------
      MODULE gas_solid_velocity
!-----------------------------------------------------
      IMPLICIT NONE
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: ug, vg, wg
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: us, vs, ws
!
      PUBLIC
      SAVE
!-----------------------------------------------------
      CONTAINS
!-----------------------------------------------------
      SUBROUTINE allocate_velocity
      USE dimensions
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE domain_mapping, ONLY: ncdom
      IMPLICIT NONE
!
      ALLOCATE(ug(ncdom), wg(ncdom))
      ALLOCATE(us(ncdom,nsolid), ws(ncdom,nsolid))
      ug = 0.0d0
      wg = 0.0d0
      us = 0.0d0
      ws = 0.0d0
      IF( job_type == JOB_TYPE_3D ) THEN
        ALLOCATE( vg(ncdom) )
        ALLOCATE( vs(ncdom,nsolid) )
        vg = 0.0d0
        vs = 0.0d0
      END IF

      RETURN
      END SUBROUTINE
!-----------------------------------------------------
      END MODULE
!-----------------------------------------------------
