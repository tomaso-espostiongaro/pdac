!----------------------------------------------------------------------
      MODULE gas_solid_temperature
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: gas_enthalpy, gas_temperature
      REAL*8, DIMENSION(:),   ALLOCATABLE :: sieg, siegn, tg
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: solid_enthalpy, solid_temperature
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: sies, siesn, ts
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE bounds_temperature
      USE dimensions
      IMPLICIT NONE
!
      ALLOCATE(gas_enthalpy(ntot), gas_temperature(ntot))
      ALLOCATE(solid_enthalpy(nsolid,ntot), solid_temperature(nsolid,ntot))
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE local_bounds_temperature
      USE dimensions
      USE grid, ONLY: ncdom, ncint
      IMPLICIT NONE
!
      ALLOCATE(sieg(ncdom), siegn(ncint), tg(ncdom))
      ALLOCATE(sies(nsolid,ncdom), siesn(nsolid,ncint))
      ALLOCATE(ts(nsolid,ncdom))
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE
