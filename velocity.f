!-----------------------------------------------------
      MODULE gas_solid_velocity
!-----------------------------------------------------
      IMPLICIT NONE
!

      REAL*8, DIMENSION(:),   ALLOCATABLE :: ug, vg, wg
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: uk, vk, wk

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
        ALLOCATE(solid_velocity_r(nsolid,ntot))
        gas_velocity_r = 0.0d0
        solid_velocity_r = 0.0d0
      ELSE
        ALLOCATE(gas_velocity_x(ntot))
        ALLOCATE(solid_velocity_x(nsolid,ntot))
        gas_velocity_x = 0.0d0
        solid_velocity_x = 0.0d0
        ALLOCATE(gas_velocity_y(ntot))
        ALLOCATE(solid_velocity_y(nsolid,ntot))
        gas_velocity_y = 0.0d0
        solid_velocity_y = 0.0d0
      END IF
      ALLOCATE(gas_velocity_z(ntot))
      ALLOCATE(solid_velocity_z(nsolid,ntot))
      gas_velocity_z = 0.0d0
      solid_velocity_z = 0.0d0
!
      RETURN
      END SUBROUTINE
!-----------------------------------------------------
      SUBROUTINE local_bounds_velocity
      USE dimensions
      USE control_flags, ONLY: job_type
      USE grid, ONLY: nijx_l
      IMPLICIT NONE
!
      ALLOCATE(ug(nijx_l), vg(nijx_l))
      ALLOCATE(uk(nsolid, nijx_l), vk(nsolid, nijx_l))
      ug = 0.0d0
      vg = 0.0d0
      uk = 0.0d0
      vk = 0.0d0
      IF( job_type == '3D' ) THEN
        ALLOCATE( wg(nijx_l) )
        ALLOCATE( wk(nsolid, nijx_l) )
        wg = 0.0d0
        wk = 0.0d0
      END IF

      RETURN
      END SUBROUTINE
!-----------------------------------------------------
      END MODULE
!-----------------------------------------------------
