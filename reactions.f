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
      SUBROUTINE hrex(ijk,hrexg,hrexs)
!
      USE dimensions
      USE time_parameters, ONLY: dt
      IMPLICIT NONE

      REAL*8 :: hrexg, hrexs
      INTEGER, INTENT(IN) :: ijk
!
      hrexs = dt * ( r1(ijk)*h1 + r2(ijk)*h2 + r3(ijk)*h3 + r4(ijk)*h4)
      hrexg = dt * (r5(ijk)*h5)
!
      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------
      SUBROUTINE rexion
!
      USE dimensions
      USE domain_decomposition, ONLY: ncint, myijk, ncdom
      USE grid, ONLY: fl_l
      USE indijk_module, ONLY: ip0_jp0_kp0_
      IMPLICIT NONE
!
      INTEGER :: ijk
      INTEGER :: imesh
!
      ALLOCATE(r1(ncdom), r2(ncdom), r3(ncdom), r4(ncdom), r5(ncdom))
!
      DO ijk = 1, ncint
          imesh = myijk( ip0_jp0_kp0_, ijk)
          IF(fl_l(ijk) == 1) THEN
            r1(ijk)=0.D0
            r2(ijk)=0.D0
            r3(ijk)=0.D0
            r4(ijk)=0.D0
            r5(ijk)=0.D0
          END IF
      END DO
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE reactions
!----------------------------------------------------------------------
