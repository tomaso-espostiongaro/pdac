!----------------------------------------------------------------------
      MODULE eos_solid
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE eosl(tk, cps, ck, siek)
!
      USE dimensions
      USE gas_constants, ONLY: tzero, hzeros
      USE specific_heat_module, ONLY: hcaps
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) :: tk
      REAL*8, INTENT(IN) ::  cps, siek
      REAL*8, INTENT(OUT) :: ck

      CALL hcaps(ck, cps, tk)
      tk = tzero + (siek-hzeros) / ck

      RETURN
      END SUBROUTINE eosl
!----------------------------------------------------------------------
      SUBROUTINE cnverts(ijk)
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
      END SUBROUTINE cnverts
!----------------------------------------------------------------------
      END MODULE eos_solid
!----------------------------------------------------------------------
