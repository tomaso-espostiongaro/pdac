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
      ALLOCATE(gas_bulk_density(ntot), gas_density(ntot))
      ALLOCATE(solid_bulk_density(nsolid,ntot))

      gas_bulk_density   = 0.D0
      gas_density        = 0.D0
      solid_bulk_density = 0.D0
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE local_bounds_density
      USE dimensions
      USE grid, ONLY: ncint, ncdom
!
      ALLOCATE(rgp(ncdom), rgpn(ncint), rog(ncint))
      ALLOCATE(rlk(ncdom,nsolid), rlkn(ncint,nsolid))
      rgp = 0.D0
      rgpn = 0.D0
      rog = 0.D0
      rlk = 0.D0
      rlkn = 0.D0
      
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE
