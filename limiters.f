!----------------------------------------------------------------------
      MODULE flux_limiters
!----------------------------------------------------------------------
      IMPLICIT NONE
!
      INTEGER, PUBLIC :: muscl
      REAL*8, PUBLIC  :: beta
      INTEGER :: lim_type
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
      REAL*8 :: beta_unlimited
      
      SELECT CASE (lim_type)

      CASE (1) !(VanLeer)

        limiter = MAX( 0.D0, (2.D0*erre/(1.D0+erre)) )

      CASE (2) !(minmod)

        limiter = MAX( 0.D0, MIN( erre, 1.D0 ))

      CASE (3) !(superbee)

        limiter = MAX( 0.D0, MIN( 2.D0*erre, 1.D0 ), MIN( erre, 2.D0 ))

      CASE (4) !(ultrabeta)

        beta_unlimited = (1.D0-beta) + beta*erre 

        limiter = MAX( 0.D0, MIN( 2.D0*erre, beta_unlimited , 2.D0 ))

      CASE (0) !(beta_unlimited)

        limiter = (1.D0-beta) + beta*erre

      CASE DEFAULT !(mantain First Order Upwind)

        limiter = 0.D0
      
      END SELECT
        
      END SUBROUTINE limiters
!-----------------------------------------------------------------------
      END MODULE flux_limiters
!-----------------------------------------------------------------------
