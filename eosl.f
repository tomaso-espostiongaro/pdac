!----------------------------------------------------------------------
      MODULE eos_solid
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE eosl(tk, cps, ck, siek, iheat,itemp)
!
      USE dimensions
      USE gas_constants, ONLY: tzero, hzeros
      USE heat_capacity, ONLY: hcaps
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iheat, itemp
      REAL*8, INTENT(INOUT) :: tk
      REAL*8, INTENT(IN) ::  cps, siek
      REAL*8 :: ck
      IF(iheat.GT.0) CALL hcaps(ck, cps, tk)
      IF(itemp.GT.0) THEN
        tk = tzero + (siek-hzeros) / ck
        CALL hcaps(ck, cps, tk)
      ENDIF
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE cnverts(ij)
!
      USE dimensions
      USE gas_solid_density, ONLY: solid_bulk_density
      USE gas_solid_temperature, ONLY: solid_enthalpy, solid_temperature
      USE gas_constants, ONLY: tzero, hzeros
      USE particles_constants, ONLY: cps
      USE reactions, ONLY: irex
      USE heat_capacity, ONLY: solid_heat_capacity, hcaps
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ij
!
      INTEGER :: k

!pdac---------------
! control next statement
!      IF(irex.EQ.0) RETURN
!pdac---------------
!
! compute heat capacity (constant volume) for particles
!
      DO k = 1, nsolid
        CALL hcaps(solid_heat_capacity(k,ij), cps(k), solid_temperature(k,ij))
        solid_enthalpy(k,ij)=(solid_temperature(k,ij)-tzero)*solid_heat_capacity(k,ij)+hzeros
      END DO
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE eos_solid
