!-----------------------------------------------------
      MODULE gas_solid_velocity
!-----------------------------------------------------
      IMPLICIT NONE
!

      REAL*8, DIMENSION(:),   ALLOCATABLE :: ug, vg
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: uk, vk

      REAL*8, DIMENSION(:),   ALLOCATABLE :: ug_g, vg_g
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: uk_g, vk_g
!
      SAVE
!-----------------------------------------------------
      CONTAINS
!-----------------------------------------------------
      SUBROUTINE bounds_velocity
      USE dimensions
      IMPLICIT NONE
!
      ALLOCATE(ug_g(ndi*ndj), vg_g(ndi*ndj))
      ALLOCATE(uk_g(ncl,ndi*ndj), vk_g(ncl,ndi*ndj))
      ug_g = 0.0d0
      vg_g = 0.0d0
      uk_g = 0.0d0
      vk_g = 0.0d0
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
      ALLOCATE(uk(ncl, nijx_l), vk(ncl, nijx_l))
      ug = 0.0d0
      vg = 0.0d0
      uk = 0.0d0
      vk = 0.0d0

      RETURN
      END SUBROUTINE
!-----------------------------------------------------
      END MODULE
!-----------------------------------------------------
