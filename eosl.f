!----------------------------------------------------------------------
      MODULE eos_solid
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE caloric_eosl(tk, cps, ck, siek, ijk)
!
      USE dimensions
      USE domain_decomposition, ONLY: meshinds
      USE gas_constants, ONLY: tzero, hzeros
      USE specific_heat_module, ONLY: hcaps
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) :: tk
      REAL*8, INTENT(IN) ::  cps, siek
      REAL*8, INTENT(OUT) :: ck
      INTEGER, INTENT(IN) :: ijk
      INTEGER :: imesh, i,j,k

      CALL hcaps( ck, cps, tk )
      tk = tzero + ( siek - hzeros ) / ck

      IF( tk < 0.0d0 ) THEN
         CALL meshinds(ijk,imesh,i,j,k)
         WRITE(6,*) 'WARNING (eosl) negative temperature'
         WRITE(6,*) 'local cell: ', ijk, i, j, k
      END IF

      RETURN
      END SUBROUTINE caloric_eosl
!----------------------------------------------------------------------
      END MODULE eos_solid
!----------------------------------------------------------------------
