!--------------------------------------------------------------------
      MODULE heat_transfer
!--------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!--------------------------------------------------------------------
      SUBROUTINE hvs(hv, rlk, rog, ep, ug, ugm, uk, ukm,
     &               vg, vgm, vk, vkm, mug, kapg, cg, k)
!
      USE dimensions
      USE particles_constants, ONLY: rl, inrl, dk 
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: k
      REAL*8, INTENT(IN) :: rlk, rog, ep
      REAL*8, INTENT(IN) :: mug, kapg, cg
      REAL*8, INTENT(IN) :: ug, ugm, uk, ukm 
      REAL*8, INTENT(IN) :: vg, vgm, vk, vkm 
      REAL*8, INTENT(OUT) :: hv
!
      REAL*8 :: reynum, pranum, eps, du, dv, vrel
      REAL*8 :: asurf
!
!pe-------------
      IF(rlk .LE. 0.D0) THEN
        hv = 0.D0
        RETURN
      ENDIF
!pe-------------
      asurf = rlk * inrl(k) * (6.D0/dk(k))
      du = 0.5D0*((ug+ugm)-(uk+ukm))
      dv = 0.5D0*((vg+vgm)-(vk+vkm))
      vrel = DSQRT(du**2+dv**2)
      reynum = dk(k) * vrel * rog / mug
      pranum = cg * mug / kapg
!
! generalized gunn's correlation for n particles
!
      eps = rlk * inrl(k)
      hv = ((2.D0 + 5.D0 * eps**2.D0)*
     $       (1.D0 + 0.7D0 * reynum**0.2D0 * pranum**(1.D0/3.D0)) +
     $       (0.13D0 + 1.2D0 * eps**2.D0) * reynum**0.7D0 *
     &        pranum**(1.D0/3.D0)) 
     &        * kapg * asurf/dk(k)
!
! gunn's correlation
!
!      hv = ((7.D0 - 10.D0 * ep + 5.D0 * ep**2.D0) *
!     $     (1.D0 + 0.7D0 * reynum**0.2 * pranum**(1.D0/3.D0)) +
!     $     (1.33D0 - 2.4D0 * ep + 1.2 * ep**2) *
!     $     reynum**0.7D0 * pranum**(1.D0/3.D0)) * kapg *
!     $     asurf/dk(k)
!
! zabrodsky's correlations
!
!      IF (ep .LE. 0.8D0) THEN
!        IF (reynum .LE. 200.D0) THEN
!          wnu = (2.D0 + 0.11D0 * reynum) * asurf
!          hv = wnu * kapg / dk(k)
!        ELSE IF (reynum .LE. 2000.D0) 
!          reynum = 4.D0 * vrel * rog/mug/asurf/dk(k)
!          hv = 0.123D0 * (reynum**0.83D0) * kapg * asurf/dk(k)
!        ELSE
!          hv = 0.61D0 * (reynum**0.67D0) * kapg * asurf/dk(k)
!        END IF
!      ELSE
!
!        IF (reynum .LE. 200.D0) THEN
!          hv = (2.D0 + 0.16D0 * reynum**0.67D0) * kapg * asurf/dk(k)
!        ELSE IF (reynum .LE. 1000.D0) THEN
!          hv = 8.2D0 * (reynum**0.6D0) * kapg * asurf/dk(k)
!        ELSE
!          hv = 1.06D0 * (reynum**0.457D0) * kapg * asurf/dk(k)
!        END IF
!      END IF
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE heat_transfer
