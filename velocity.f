!-----------------------------------------------------
      MODULE gas_solid_velocity
!-----------------------------------------------------
      IMPLICIT NONE
!

      REAL*8, DIMENSION(:),   ALLOCATABLE :: ug, vg, wg
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: us, vs, ws

      REAL*8, DIMENSION(:),   ALLOCATABLE :: gas_velocity_r, gas_velocity_z, gas_velocity_x, gas_velocity_y
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: solid_velocity_r, solid_velocity_z, solid_velocity_x, solid_velocity_y
!
      SAVE
!-----------------------------------------------------
      CONTAINS
!-----------------------------------------------------
      SUBROUTINE bounds_velocity
      USE dimensions
      USE control_flags, ONLY: job_type
      IMPLICIT NONE
!
      IF( job_type == '2D' ) THEN
        ALLOCATE(gas_velocity_r(ntot))
        ALLOCATE(solid_velocity_r(ntot,nsolid))
        gas_velocity_r = 0.0d0
        solid_velocity_r = 0.0d0
      ELSE
        ALLOCATE(gas_velocity_x(ntot))
        ALLOCATE(solid_velocity_x(ntot,nsolid))
        gas_velocity_x = 0.0d0
        solid_velocity_x = 0.0d0
        ALLOCATE(gas_velocity_y(ntot))
        ALLOCATE(solid_velocity_y(ntot,nsolid))
        gas_velocity_y = 0.0d0
        solid_velocity_y = 0.0d0
      END IF
      ALLOCATE(gas_velocity_z(ntot))
      ALLOCATE(solid_velocity_z(ntot,nsolid))
      gas_velocity_z = 0.0d0
      solid_velocity_z = 0.0d0
!
      RETURN
      END SUBROUTINE
!-----------------------------------------------------
      SUBROUTINE local_bounds_velocity
      USE dimensions
      USE control_flags, ONLY: job_type
      USE grid, ONLY: ncdom
      IMPLICIT NONE
!
      ALLOCATE(ug(ncdom), wg(ncdom))
      ALLOCATE(us(ncdom,nsolid), ws(ncdom,nsolid))
      ug = 0.0d0
      wg = 0.0d0
      us = 0.0d0
      ws = 0.0d0
      IF( job_type == '3D' ) THEN
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
