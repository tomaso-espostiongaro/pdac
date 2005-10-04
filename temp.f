!----------------------------------------------------------------------
      MODULE gas_solid_temperature
!----------------------------------------------------------------------
      IMPLICIT NONE
!
! ... gas specific internal enthalpy, gas temperature

      REAL*8, DIMENSION(:),   ALLOCATABLE :: sieg, siegn, tg
!
! ... solid specific internal enthalpy, solid temperature

      REAL*8, DIMENSION(:,:), ALLOCATABLE :: sies, siesn, ts
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_temperature
      USE dimensions
      USE domain_mapping, ONLY: ncdom, ncint
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
