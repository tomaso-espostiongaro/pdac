!----------------------------------------------------------------------
      MODULE gas_solid_temperature
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: gas_enthalpy, gas_temperature
      REAL*8, DIMENSION(:),   ALLOCATABLE :: sieg, siegn, tg
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: solid_enthalpy, solid_temperature
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: siek, siekn, tk
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE bounds_temperature
      USE dimensions
      IMPLICIT NONE
!
      ALLOCATE(gas_enthalpy(nr*nz), gas_temperature(nr*nz))
      ALLOCATE(solid_enthalpy(nsolid,nr*nz), solid_temperature(nsolid,nr*nz))
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE local_bounds_temperature
      USE dimensions
      USE grid, ONLY: nijx_l, nij_l
      IMPLICIT NONE
!
      ALLOCATE(sieg(nijx_l), siegn(nijx_l), tg(nijx_l))
      ALLOCATE(siek(nsolid,nijx_l), siekn(nsolid,nijx_l))
      ALLOCATE(tk(nsolid,nijx_l))
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE
