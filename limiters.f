!----------------------------------------------------------------------
      MODULE flux_limiters
!----------------------------------------------------------------------
      IMPLICIT NONE
!
      INTEGER, PUBLIC :: muscl
      REAL*8, PUBLIC  :: beta

      LOGICAL, PRIVATE:: vanleer, minmod, superbee, ultrabeta
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE limiters(limiter,erre)
! ... this routine computes the flux limiter accordingly to the
! ... different methods reviewed in Sweby'84' and Leonard'90'
! ... ('Ultra-Sharp' technique) papers
      IMPLICIT NONE

      REAL*8, INTENT(OUT) :: limiter
      REAL*8, INTENT(IN)  :: erre
      
      vanleer   = .FALSE.
      minmod    = .FALSE.
      superbee  = .FALSE.
      ultrabeta = .TRUE.

      IF (vanleer) THEN
        limiter = ( DABS(erre) + erre ) / ( 1.D0 + DABS(erre) )
      ELSE IF (minmod) THEN
        limiter = MAX( 0.D0, MIN( erre, 1.D0 ))
      ELSE IF (superbee) THEN
        limiter = MAX( 0.D0, MIN( 2.D0*erre, 1.D0 ), MIN( erre, 2.D0 ))
      ELSE IF (ultrabeta) THEN
        limiter = MAX( 0.D0, MIN( 2.D0*erre, ((1.D0-beta) + beta*erre), 2.D0 ))
      END IF
        
      END SUBROUTINE limiters
!-----------------------------------------------------------------------
      END MODULE flux_limiters
!-----------------------------------------------------------------------
