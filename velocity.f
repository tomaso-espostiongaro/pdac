!-----------------------------------------------------
      MODULE gas_solid_velocity
!-----------------------------------------------------
      IMPLICIT NONE
!

      REAL*8, DIMENSION(:),   ALLOCATABLE :: ug, vg
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: uk, vk

      REAL*8, DIMENSION(:),   ALLOCATABLE :: gas_velocity_r, gas_velocity_z
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: solid_velocity_r, solid_velocity_z
!
      SAVE
!-----------------------------------------------------
      CONTAINS
!-----------------------------------------------------
      SUBROUTINE bounds_velocity
      USE dimensions
      IMPLICIT NONE
!
      ALLOCATE(gas_velocity_r(nr*nz), gas_velocity_z(nr*nz))
      ALLOCATE(solid_velocity_r(nsolid,nr*nz), solid_velocity_z(nsolid,nr*nz))
      gas_velocity_r = 0.0d0
      gas_velocity_z = 0.0d0
      solid_velocity_r = 0.0d0
      solid_velocity_z = 0.0d0
!
      RETURN
      END SUBROUTINE
!-----------------------------------------------------
      SUBROUTINE local_bounds_velocity
      USE dimensions
      USE grid, ONLY: nijx_l
      IMPLICIT NONE
!
      ALLOCATE(ug(nijx_l), vg(nijx_l))
      ALLOCATE(uk(nsolid, nijx_l), vk(nsolid, nijx_l))
      ug = 0.0d0
      vg = 0.0d0
      uk = 0.0d0
      vk = 0.0d0

      RETURN
      END SUBROUTINE
!-----------------------------------------------------
      END MODULE
!-----------------------------------------------------
