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
      USE specific_heat_module, ONLY: hcaps
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
      SUBROUTINE cnverts_local(ijk)
!
      USE dimensions
      USE gas_solid_temperature, ONLY: sies, ts
      USE gas_constants, ONLY: tzero, hzeros
      USE particles_constants, ONLY: cps
      USE specific_heat_module, ONLY: ck, hcaps
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ijk
      INTEGER :: is
!
! compute heat capacity (constant volume) for particles
!
      DO is = 1, nsolid
        CALL hcaps(ck(is,ijk), cps(is), ts(ijk,is))
        sies(ijk,is) = ( ts(ijk,is) - tzero ) * ck(is,ijk) + hzeros
      END DO
!
      RETURN
      END SUBROUTINE cnverts_local


!----------------------------------------------------------------------
      END MODULE eos_solid
!----------------------------------------------------------------------
