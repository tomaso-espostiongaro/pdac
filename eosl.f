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
      USE specific_heat, ONLY: hcaps
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iheat, itemp
      REAL*8, INTENT(INOUT) :: tk
      REAL*8, INTENT(IN) ::  cps, siek
      REAL*8, INTENT(OUT) :: ck

      IF(iheat > 0) CALL hcaps(ck, cps, tk)
      IF(itemp > 0) THEN
        tk = tzero + (siek-hzeros) / ck
        CALL hcaps(ck, cps, tk)
      ENDIF

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE cnverts(imesh)
!
      USE dimensions
      USE gas_solid_temperature, ONLY: solid_enthalpy, solid_temperature
      USE gas_constants, ONLY: tzero, hzeros
      USE particles_constants, ONLY: cps
      USE specific_heat, ONLY: solid_specific_heat, hcaps
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: imesh
      INTEGER :: is
!
! compute heat capacity (constant volume) for particles
!
      DO is = 1, nsolid
        CALL hcaps(solid_specific_heat(is,imesh), cps(is), &
                   solid_temperature(imesh,is))
        solid_enthalpy(imesh,is) = ( solid_temperature(imesh,is) - tzero ) * &
                                    solid_specific_heat(is,imesh) + hzeros
      END DO
!
      RETURN
      END SUBROUTINE cnverts
!----------------------------------------------------------------------
      END MODULE eos_solid
!----------------------------------------------------------------------
