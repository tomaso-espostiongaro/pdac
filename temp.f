!----------------------------------------------------------------------
      MODULE gas_solid_temperature
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: sieg_g, siegn_g, tg_g
      REAL*8, DIMENSION(:),   ALLOCATABLE :: sieg, siegn, tg
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: siek_g, siekn_g, tk_g
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: siek, siekn, tk
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE bounds_temperature
      USE dimensions
      IMPLICIT NONE
!
      ALLOCATE(sieg_g(ndi*ndj), siegn_g(ndi*ndj), tg_g(ndi*ndj))
      ALLOCATE(siek_g(ncl,ndi*ndj), siekn_g(ncl,ndi*ndj), tk_g(ncl,ndi*ndj))
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE local_bounds_temperature
      USE dimensions
      USE grid, ONLY: nijx_l, nij_l
      IMPLICIT NONE
!
      ALLOCATE(sieg(nijx_l), siegn(nijx_l), tg(nijx_l))
      ALLOCATE(siek(ncl,nijx_l), siekn(ncl,nijx_l))
      ALLOCATE(tk(ncl,nijx_l))
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE
