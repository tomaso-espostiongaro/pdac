!----------------------------------------------------------------------
      MODULE gas_solid_density
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: rgp_g, rgpn_g, rog_g
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rlk_g, rlkn_g

      REAL*8, DIMENSION(:), ALLOCATABLE :: rgp, rgpn, rog
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rlk, rlkn
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE bounds_density
      USE dimensions
!
      ALLOCATE(rgp_g(ndi*ndj), rgpn_g(ndi*ndj), rog_g(ndi*ndj))
      ALLOCATE(rlk_g(ncl,ndi*ndj), rlkn_g(ncl,ndi*ndj))

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE local_bounds_density
      USE dimensions
      USE grid, ONLY: nij_l, nijx_l
!
      ALLOCATE(rgp(nijx_l), rgpn(nij_l), rog(nij_l))
      ALLOCATE(rlk(ncl,nijx_l), rlkn(ncl,nij_l))
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE
