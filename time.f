!---------------------------------------------------
      MODULE time_parameters
!---------------------------------------------------
      IMPLICIT NONE
      SAVE

      INTEGER :: itd, rungekut, sweep, time_integration
      INTEGER :: ndump, nprint
      REAL*8 :: time, tdump, tpr, tstop, dt, timestart
      REAL*8 :: tau, tau1, tau2
      INTEGER :: ift
      REAL*8 :: alpha, alphagrav
      LOGICAL :: adapt_dt
!---------------------------------------------------
      END MODULE
!---------------------------------------------------
