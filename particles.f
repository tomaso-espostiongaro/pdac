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
        REAL*8, DIMENSION(:),   ALLOCATABLE :: plim           ! limit for momentum matrix
        REAL*8, DIMENSION(:,:), ALLOCATABLE :: dkf            ! p-p interaction coefficient
        REAL*8, DIMENSION(:,:), ALLOCATABLE :: philim
        REAL*8, DIMENSION(:,:), ALLOCATABLE :: epsl
        REAL*8, DIMENSION(:,:), ALLOCATABLE :: epsu
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_part_constants
      USE dimensions
!
      ALLOCATE(dk(nsolid), dkm1(nsolid), rl(nsolid), rlm1(nsolid),      &
               phis(nsolid), phism1(nsolid), cmus(nsolid), cps(nsolid), &
               phi(nsolid), kap(nsolid), inrl(nsolid), plim(nsolid))
      ALLOCATE(dkf(nsolid,nsolid), philim(nsolid,nsolid),               &
               epsl(nsolid,nsolid), epsu(nsolid,nsolid))
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE particles_constants_set
!
      USE dimensions
      IMPLICIT NONE
!
      INTEGER :: k2, k1, k, kk
      REAL*8  :: a, fact
!
! ... set the maximum packing volume fraction for particles
! ... (here assumed to be the same for all phases)
!
      phi(:) = 0.63
!
! ... Syamlal's particle-particle interaction coefficients
!
      DO k=1,nsolid
        DO kk=1,nsolid
          IF(dk(k) >= dk(kk)) THEN
            k1=k
            k2=kk
          ELSE
            k1=kk
            k2=k
          ENDIF
          a = DSQRT(dk(k2)/dk(k1))
          
          philim(k,kk) = phi(k1) / ( phi(k1) + ( 1.D0 - phi(k1)) * phi(k2) )

          epsl(k,kk) = (phi(k1)-phi(k2))+((1.D0-a)*(1.D0-phi(k1))*phi(k2)) 

          fact = ( phi(k1) + ( 1.D0 - phi(k2) ) * phi(k1) ) / phi(k1)

          epsl(k,kk) = epsl(k,kk) * fact

          epsu(k,kk) = ( 1.D0 - a ) * ( phi(k1) + ( 1.D0 - phi(k1) ) * phi(k2) )
!
          dkf(k,kk)=(dk(k)+dk(kk))**2/(rl(k)*dk(k)**3+rl(kk)*dk(kk)**3)
          dkf(kk,k)=dkf(k,kk)
!
        END DO
      END DO
!
      RETURN
      END SUBROUTINE particles_constants_set
!----------------------------------------------------------------------
      END MODULE particles_constants
!----------------------------------------------------------------------
