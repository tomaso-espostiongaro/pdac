!----------------------------------------------------------------------
      MODULE eos_solid
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE caloric_eosl(tk, cps, ck, siek)
!
      USE dimensions
      USE gas_constants, ONLY: tzero, hzeros
      USE specific_heat_module, ONLY: hcaps
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) :: tk
      REAL*8, INTENT(IN) ::  cps, siek
      REAL*8, INTENT(OUT) :: ck

      CALL hcaps( ck, cps, tk )
      tk = tzero + ( siek - hzeros ) / ck

      IF( tk < 0.0d0 ) THEN
         WRITE(6,*) 'WARNING (eosl) negative temperature'
         tk = 0.0d0
      END IF

      RETURN
      END SUBROUTINE caloric_eosl
!----------------------------------------------------------------------
      END MODULE eos_solid
!----------------------------------------------------------------------
