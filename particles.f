!----------------------------------------------------------------------
      MODULE particles_constants
!----------------------------------------------------------------------
        USE dimensions
        IMPLICIT NONE
        SAVE

        REAL*8, DIMENSION(:),   ALLOCATABLE :: dk, dkm1       ! particle diameters
        REAL*8, DIMENSION(:),   ALLOCATABLE :: rl, rlm1       ! specIFic density
        REAL*8, DIMENSION(:),   ALLOCATABLE :: inrl           ! inverse of specIFic density
        REAL*8, DIMENSION(:),   ALLOCATABLE :: phis, phism1   ! sphericity
        REAL*8, DIMENSION(:),   ALLOCATABLE :: cmus
        REAL*8, DIMENSION(:),   ALLOCATABLE :: cps            ! thermal capacity
        REAL*8, DIMENSION(:),   ALLOCATABLE :: phi
        REAL*8, DIMENSION(:),   ALLOCATABLE :: kap            ! thermal conductivity
        REAL*8, DIMENSION(:,:), ALLOCATABLE :: dkf            ! p-p interaction coefficient
        REAL*8, DIMENSION(:,:), ALLOCATABLE :: philim
        REAL*8, DIMENSION(:,:), ALLOCATABLE :: epsl
        REAL*8, DIMENSION(:,:), ALLOCATABLE :: epsu
        INTEGER :: nsolid
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE bounds_part_constants
      USE dimensions
!
      ALLOCATE(dk(ncl), dkm1(ncl), rl(ncl), rlm1(ncl), phis(ncl), phism1(ncl), &
               cmus(ncl), cps(ncl), phi(ncl), kap(ncl), inrl(ncl))
      ALLOCATE(dkf(ncl,ncl), philim(ncl,ncl), epsl(ncl,ncl), epsu(ncl,ncl))
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE particles_constants_set(nsolid)
!
      USE dimensions
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: nsolid
!
      INTEGER :: k2, k1, k, kk
      REAL*8 :: dratx

!
! Syamlal's particle-particle interaction coefficients
!
      DO k=1,nsolid
        phi(k)=0.63
      END DO
!
      DO k=1,nsolid
      DO kk=1,nsolid
        IF(dk(k).GE.dk(kk)) THEN
          k1=k
          k2=kk
        ELSE
          k1=kk
          k2=k
        ENDIF
        dratx=DSQRT(dk(k2)/dk(k1))
        philim(k,kk)=phi(k1)/(phi(k1)+(1.D0-phi(k1))*phi(k2))
        epsl(k,kk)=(phi(k1)-phi(k2)+(1.D0-dratx)*(1.D0-phi(k1))*phi(k2))  &
                   *(phi(k1)+(1.D0-phi(k2))*phi(k1))
        epsu(k,kk)=(1.D0-dratx)*(phi(k1)+(1.D0-phi(k1))*phi(k2))
      END DO
      END DO
!
      DO k=1,nsolid
        DO kk=1,k
          dkf(k,kk)=(dk(k)+dk(kk))**2/(rl(k)*dk(k)**3+rl(kk)*dk(kk)**3)
          dkf(kk,k)=dkf(k,kk)
        END DO
      END DO
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE particles_constants

