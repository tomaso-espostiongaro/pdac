!----------------------------------------------------------------------
      MODULE reactions
!----------------------------------------------------------------------
! ... DUMMY REACTION MODULE ...
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
      USE grid, ONLY: nij_l, myij, nijx_l
      USE dimensions
      IMPLICIT NONE
!
      INTEGER :: ij
      INTEGER :: imesh
!
      ALLOCATE(r1(nijx_l), r2(nijx_l), r3(nijx_l), r4(nijx_l), r5(nijx_l))
!
      DO ij = 1, nij_l
          imesh = myij(0, 0, ij)
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
