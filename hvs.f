!--------------------------------------------------------------------
      MODULE heat_transfer
!--------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!--------------------------------------------------------------------
      SUBROUTINE hvs(hv, rlk, rog, ep, du, dv, dw, mug, kapg, cg, k)
!
      USE dimensions
      USE particles_constants, ONLY: rl, inrl, dk 
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: k
      REAL*8, INTENT(IN) :: rlk, rog, ep
      REAL*8, INTENT(IN) :: mug, kapg, cg
      REAL*8, INTENT(IN) :: du, dv, dw
      REAL*8, INTENT(OUT) :: hv
!
      REAL*8 :: reynum, pranum, eps, eps2, vrel
      REAL*8 :: asurf
      REAL*8 :: nusselt
!
!pe-------------
      IF(rlk <= 0.D0) THEN
        hv = 0.D0
        RETURN
      ENDIF
!pe-------------
      eps  =  rlk * inrl(k)
      eps2 = (rlk * inrl(k))**2
      asurf = 6.D0 * kapg * eps / (dk(k)**2)
      vrel = DSQRT(du**2+dv**2+dw**2)
      reynum = dk(k) * vrel * rog / mug
      pranum = cg * mug / kapg
!
! generalized Gunn's correlation for n particles
!
      nusselt = ((2.D0 + 5.D0 * eps2)*                                   &
                (1.D0 + 0.7D0 * reynum**0.2D0 * pranum**(1.D0/3.D0)) +   &
                (0.13D0 + 1.2D0 * eps2) * reynum**0.7D0 *                &
                pranum**(1.D0/3.D0))
!
! Gunn's correlation
!
!      nusselt = ((7.D0 - 10.D0 * ep + 5.D0 * ep**2) *                &
!                (1.D0 + 0.7D0 * reynum**0.2 * pranum**(1.D0/3.D0)) +    &
!                (1.33D0 - 2.4D0 * ep + 1.2 * ep**2) *                   &
!                reynum**0.7D0 * pranum**(1.D0/3.D0))
!
! Zabrodsky's correlations
!
!      IF (ep <= 0.8D0) THEN
!        IF (reynum <= 200.D0) THEN
!          nusselt = (2.D0 + 0.11D0 * reynum)
!        ELSE IF (reynum <= 2000.D0) THEN
!          nusselt = 0.123D0 * (2.D0 * reynum / (3.D0 * eps))**0.83D0
!        ELSE IF (reynum > 2000.D0) THEN
!          nusselt = 0.61D0 * reynum**0.67D0
!        END IF
!      ELSE IF (ep > 0.8D0) THEN
!        IF (reynum <= 200.D0) THEN
!          nusselt = (2.D0 + 0.16D0 * reynum**0.67D0)
!        ELSE IF (reynum <= 1000.D0) THEN
!          nusselt = 8.2D0 * (reynum**0.6D0)
!        ELSE IF (reynum > 1000.D0) THEN
!          nusselt = 1.06D0 * (reynum**0.457D0)
!        END IF
!      END IF
!
      hv = nusselt * asurf

      RETURN
      END SUBROUTINE hvs
!----------------------------------------------------------------------
      END MODULE heat_transfer
!----------------------------------------------------------------------
