!----------------------------------------------------------------------
      MODULE eos_solid
!----------------------------------------------------------------------
      USE io_files, ONLY: testunit
      IMPLICIT NONE
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE caloric_eosl(tk, cps, ck, siek, ijk, info)
!
      USE control_flags, ONLY: lpr
      USE dimensions
      USE domain_mapping, ONLY: meshinds
      USE gas_constants, ONLY: tzero, hzeros
      USE parallel, ONLY: mpime
      USE specific_heat_module, ONLY: hcaps
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) :: tk
      REAL*8, INTENT(IN) ::  cps, siek
      REAL*8, INTENT(OUT) :: ck
      INTEGER, INTENT(IN) :: ijk
      INTEGER :: imesh, i,j,k
      INTEGER, INTENT(INOUT) :: info

      CALL hcaps( ck, cps, tk )
      tk = tzero + ( siek - hzeros ) / ck

      IF( tk < 0.0d0) THEN
        info = info + 1
        IF(lpr>=1) THEN
           CALL meshinds(ijk,imesh,i,j,k)
           WRITE(testunit,*) 'WARNING from proc: ', mpime
           WRITE(testunit,*) 'negative temperature in eosl'
           WRITE(testunit,*) 'local cell: ', ijk, i, j, k
        END IF
      END IF

      RETURN
      END SUBROUTINE caloric_eosl
!----------------------------------------------------------------------
      END MODULE eos_solid
!----------------------------------------------------------------------
