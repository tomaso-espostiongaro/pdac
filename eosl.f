!----------------------------------------------------------------------
      MODULE eos_solid
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE eosl(tk, cps, ck, siek, nc,nt)
!
      USE dimensions
      USE gas_constants, ONLY: tzero, hzeros
      USE th_capacity, ONLY: hcaps
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nc, nt
      REAL*8, INTENT(INOUT) :: tk
      REAL*8, INTENT(IN) ::  cps, siek
      REAL*8 :: ck
      IF(nc.GT.0) CALL hcaps(ck, cps, tk)
      IF(nt.GT.0) THEN
        tk = tzero + (siek-hzeros) / ck
        CALL hcaps(ck, cps, tk)
      ENDIF
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE cnverts(ij)
!
      USE dimensions
      USE gas_solid_density, ONLY: rlk_g, rlkn_g
      USE gas_solid_temperature, ONLY: siek_g, siekn_g, tk_g
      USE gas_constants, ONLY: tzero, hzeros
      USE particles_constants, ONLY: nsolid, cps
      USE reactions, ONLY: irex
      USE th_capacity, ONLY: ck_g, hcaps
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ij
!
      REAL*8 :: c1, c2
      INTEGER :: k

! ... densita delle particelle per cella
      DO k = 1, nsolid
        rlkn_g(k,ij) = rlk_g(k,ij)
      END DO

!pdac---------------
! control next statement
!      IF(irex.EQ.0) RETURN
!pdac---------------
!
! compute heat capacity (constant volume) for particles
!
      DO k = 1, nsolid
        CALL hcaps(ck_g(k,ij), cps(k), tk_g(k,ij))
        siek_g(k,ij)=(tk_g(k,ij)-tzero)*ck_g(k,ij)+hzeros
        siekn_g(k,ij)=siek_g(k,ij)
      END DO
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE eos_solid
