!----------------------------------------------------------------------
      MODULE gas_solid_density
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: gas_density, gas_bulk_density
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: solid_bulk_density

      REAL*8, DIMENSION(:), ALLOCATABLE :: rgp, rgpn, rog
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rlk, rlkn
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE bounds_density
      USE dimensions
!
      ALLOCATE(gas_bulk_density(ndi*ndj), gas_density(ndi*ndj))
      ALLOCATE(solid_bulk_density(ncl,ndi*ndj))

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
