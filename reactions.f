!----------------------------------------------------------------------
      MODULE reactions
!----------------------------------------------------------------------
! ... DUMMY REACTION MODULE ...
!
!
      IMPLICIT NONE
      INTEGER :: irex
      REAL*8, DIMENSION(:), ALLOCATABLE :: r1, r2, r3, r4, r5
      REAL*8 :: h1, h2, h3, h4, h5 
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE hrex(ij,hrexg,hrexs)
!
      USE dimensions
      USE time_parameters, ONLY: dt
      IMPLICIT NONE

      REAL*8 :: hrexg, hrexs
      INTEGER, INTENT(IN) :: ij
!
      hrexs = dt * ( r1(ij)*h1 + r2(ij)*h2 + r3(ij)*h3 + r4(ij)*h4)
      hrexg = dt * (r5(ij)*h5)
!
      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------
      SUBROUTINE rexion
!
      USE grid, ONLY: fl_l
      USE grid, ONLY: ncint, myijk, ncdom
      USE dimensions
      USE indijk_module, ONLY: ip0_jp0_kp0_
      IMPLICIT NONE
!
      INTEGER :: ij
      INTEGER :: imesh
!
      ALLOCATE(r1(ncdom), r2(ncdom), r3(ncdom), r4(ncdom), r5(ncdom))
!
      DO ij = 1, ncint
          imesh = myijk( ip0_jp0_kp0_, ij)
          IF(fl_l(ij).EQ.1) THEN
            r1(ij)=0.D0
            r2(ij)=0.D0
            r3(ij)=0.D0
            r4(ij)=0.D0
            r5(ij)=0.D0
          END IF
      END DO
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE reactions
