!------------------------------------------------------------------------
      MODULE atmosphere
!------------------------------------------------------------------------
        IMPLICIT NONE
        SAVE
!
! ... Initial atmospheric variables
!
        REAL*8 :: u0, v0, w0, p0, temp0, us0, vs0, ws0
        REAL*8 :: ep0, epsmx0
        REAL*8 :: gravx, gravz
        LOGICAL :: stratification
!------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------
      SUBROUTINE controlatm

      stratification = .TRUE.
      IF (gravz == 0.D0) THEN
         stratification = .FALSE.
      ELSE IF (gravz /= -9.81D0) THEN
         WRITE(*,*) 'WARNING!! control atmospheric stratification'
      END IF
!
      IF (stratification) THEN
        IF (temp0 /= 288.15D0) WRITE(*,*) 'WARNING! control atmospheric &
                                           temperature profile'
        IF (p0 /= 1.01325D5) WRITE(*,*)   'WARNING! control atmospheric &
                                           pressure profile'
      END IF	
!
      END SUBROUTINE controlatm
!------------------------------------------------------------------------
      SUBROUTINE atm( za, pa, ta )
!------------------------------------------------------------------------
! ... This routine computes atmospheric stratification accordingly with
! ... teh description of a standard, quiet atmosphere
!
      USE dimensions
      USE gas_constants, ONLY: gmw,rgas
!
      IMPLICIT NONE

      REAL*8, INTENT(IN) :: za 
      REAL*8, INTENT(OUT) :: pa, ta
      REAL*8 :: p1, t1, p2, t2, p3, t3, p4, t4, p5, t5, p6, t6
!
      REAL*8 :: erreair, gradt, cost
!
      erreair=rgas/gmw(6)
!
!pe------------------------------
! zone 0
      IF(za <= 0.D0 .OR. .NOT.stratification) THEN
        ta=temp0
        pa=p0
!pe------------------------------
!
! zone 1: constant (negative) temperature gradient (Troposphere)
!
      ELSE IF(za <= 1.1D4) THEN
        gradt = - 6.5D-3
        cost = gravz/(erreair*gradt)
        ta    = temp0 + gradt*za
        pa    = p0 * (ta/temp0)**cost
!
! zone 2: constant temperature (Tropopause)
!
      ELSE IF(za <= 2.0D4) THEN
	p1    = 22620.45D0
	t1    = 216.65D0
        ta    = 216.65D0
        cost  = gravz/(erreair*ta)
        pa    = p1*DEXP(cost*(za-1.1D4))
!
! zone 3: constant (positive) temperature gradient (Stratosphere)
!
      ELSE IF(za <= 3.2D4) THEN
        p2    = 5477.79D0
	t2    = 216.65D0
        gradt = 1.D-3
        cost  = gravz/(erreair*gradt)
        ta    = t2 + gradt*(za - 2.0D4)
        pa    = p2*(ta/t2)**cost
!
! zone 4: constant (positive) temperature gradient (Stratosphere)
!
      ELSE IF(za <= 4.7D4) THEN
        p3    = 869.19D0
	t3    = 228.65D0
	gradt = 2.8D-5
        cost  = gravz/(erreair*gradt)
        ta    = t3 + gradt*(za - 3.2D4)
        pa    = p3*(ta/t3)**cost
!
! zone 5: constant temperature (Ozone layer)
!
      ELSE IF(za <= 5.1D4) THEN
	p4    = 111.19D0
	t4    = 270.65D0
        ta    = 270.65D0
        cost  = gravz/(erreair*ta)
        pa    = p4*DEXP(cost*(za-4.7D4))
!
! zone 6: constant (negative) temperature gradient (Mesosphere)
!
        ELSE IF(za <= 7.1D4) THEN
	p5    = 67.24D0
	t5    = 270.65D0
	gradt = - 2.8D-3
        cost  = gravz/(erreair*gradt)
        ta    = t5 + gradt*(za - 5.1D4)
        pa    = p5*(ta/t5)**cost
!
! zone 7: constant(negative) temperature gradient (Mesosphere)
!
      ELSE IF(za <= 8.0D4) THEN
	p6    = 3.97D0
	t6    = 214.65D0
	gradt = - 2.0D-3
        cost = gravz/(erreair*gradt)
        ta    = t6 + gradt*(za - 7.1D4)
        pa    = p6*(ta/t6)**cost
!
      ELSE
        CALL error('atm','altitude out of range',1)
      ENDIF
!
      RETURN
      END SUBROUTINE atm
!----------------------------------------------------------------------
      END MODULE atmosphere
!------------------------------------------------------------------------
