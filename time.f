!---------------------------------------------------
      MODULE time_parameters
!---------------------------------------------------
      IMPLICIT NONE
      SAVE

      INTEGER :: itd, rungekut, sweep
      INTEGER :: ndump, nprint
      REAL*8 :: time, tdump, tpr, tstop, dt, timestart
      REAL*8 :: tau1, tau2
!---------------------------------------------------
      END MODULE
!---------------------------------------------------
