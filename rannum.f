!----------------------------------------------------------------------
! ... Minimal Standard Portable Random Number Generator
! ... (Numerical Recipes in Fortran, 1992)
!----------------------------------------------------------------------
      REAL*8 FUNCTION ran0(idum)
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: idum
      INTEGER :: ia, im, iq, ir, mask
      REAL*8  :: am
      PARAMETER( ia = 16807, im = 2147483647, am = 1.D0/im,  &
                 iq = 127773, ir = 2836, mask = 123459876 )
      INTEGER :: k
      
      idum = IEOR(idum, mask)
      k = idum / iq
      idum = ia * ( idum - k*iq ) - ir*k
      IF (idum < 0) idum = idum + im
      ran0 = am * idum
      idum = IEOR(idum, mask)

      RETURN
      END FUNCTION ran0
!----------------------------------------------------------------------
