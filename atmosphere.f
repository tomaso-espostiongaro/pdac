!------------------------------------------------------------------------
      MODULE atmosphere
!------------------------------------------------------------------------
        IMPLICIT NONE
        SAVE
!
! ... Initial atmospheric variables
!
        REAL*8 :: u0, v0, p0, temp0, uk0, vk0
        REAL*8 :: ep0, epsmx0
        REAL*8 :: gravx, gravz
        LOGICAL :: stratification
!------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------
      SUBROUTINE atm(za,pa,ta)
!------------------------------------------------------------------------
! ... This routine computes atmospheric stratification accordingly with
! ... teh description of a standard, quiet atmosphere
!
      USE dimensions
      USE gas_constants, ONLY: gmw,rgas
!
      IMPLICIT NONE

      REAL*8 :: za, pa, ta
!
      REAL*8 :: erreair
      REAL*8 :: cost1, cost2, cost3, cost4, cost5, cost6, cost7
!
      erreair=rgas/gmw(6)
!
      stratification = .FALSE.
!
! sostituire le costanti numeriche con dei PARAMETER !!
!
!pe------------------------------
! zone 0
      IF(za.LE.0D0 .OR. .NOT.stratification) THEN
        ta=temp0
        pa=p0
!pe------------------------------
! zone 1
      ELSE IF(za.LE.1.1D6) THEN
        cost1=-gravz/(erreair*6.5D-5)
        ta=temp0-6.5D-5*za
        pa=p0*((temp0-za*6.5D-5)/temp0)**cost1
! zone 2
      ELSE IF(za.LE.2.D6) THEN
        cost2=-gravz/(erreair*216.65D0)
        ta=216.65D0
        pa=226550.5695D0/DEXP((za-1.1D6)*cost2)
! zone 3
      ELSE IF(za.LE.3.2D6) THEN
        cost3=-gravz/(erreair*1.D-5)
        ta=196.65D0+za*1.D-5
        pa=54857.2227D0*(216.65D0/(196.65D0+za*1.D-5))**cost3
! zone 4
      ELSE IF(za.LE.4.7D6) THEN
        cost4=-gravz/(erreair*2.8D-5)
        ta=139.05D0+za*2.8D-5
        pa=8708.2208D0*(228.65D0/(139.05D0+za*2.8D-5))**cost4
! zone 5
      ELSE IF(za.LE.5.1D6) THEN
        cost5=-gravz/(erreair*270.65D0)
        ta=270.65D0
        pa=1114.1968D0/dEXP((za-4.7D6)*cost5)
! zone 6
      ELSE IF(za.LE.7.1D6) THEN
        cost6=-gravz/(erreair*2.8D-5)
        ta=413.45D0-2.8D-5*za
        pa=672.7173D0*((413.45D0-za*2.8D-5)/270.65D0)**cost6
! zone 7
      ELSE IF(za.LE.8.0D6) THEN
        cost7=-gravz/(erreair*2.0D-5)
        ta=356.65D0-2.0D-5*za
        pa=39.8372D0*((356.65D0-za*2.0D-5)/214.65D0)**cost7
      ELSE
        WRITE(8,*) 'altitude out of range'
        STOP
      ENDIF
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE
!------------------------------------------------------------------------
