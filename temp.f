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
      SUBROUTINE allocate_temperature
      USE dimensions
      USE grid, ONLY: ncdom, ncint
      IMPLICIT NONE
!
      ALLOCATE(sieg(ncdom), siegn(ncint), tg(ncdom))
      ALLOCATE(sies(ncdom,nsolid), siesn(ncint,nsolid))
      ALLOCATE(ts(ncdom,nsolid))
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE
!----------------------------------------------------------------------
