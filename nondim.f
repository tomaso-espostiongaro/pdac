!------------------------------------------------------------------------
      MODULE nondim_numbers
!------------------------------------------------------------------------
      USE dimensions
      USE domain_decomposition, ONLY: ncint
      IMPLICIT NONE
      SAVE
!
      REAL*8, ALLOCATABLE :: rich(:)        ! Richardson number
      REAL*8, ALLOCATABLE :: mut2mu(:)      ! Turbulent/Molecular viscosity
      REAL*8, ALLOCATABLE :: gas2pp(:,:)    ! Gas-Particle/Particle-Particle
      REAL*8, ALLOCATABLE :: drag2grav(:,:) ! Gravity/Gas-Particle
!------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------
      SUBROUTINE drag_ratio(appu,appw,kpgv,rlk,ijk)
! ... compute the ratio between gas-particle drag and the sum of particle-particle 
! ... drag coefficients for each solid phase
!
      USE time_parameters, ONLY: time, dt
      IMPLICIT NONE

      REAL*8,  INTENT(IN) :: appu(:), appw(:)
      REAL*8,  INTENT(IN) :: kpgv(:), rlk(:)
      INTEGER, INTENT(IN) :: ijk
      INTEGER :: is, n, ks, ks1
      REAL*8 :: gasdrag, ppdrag, ppx, ppy
!
      DO is = 1, nsolid
        n  = is + 1
        ks = (n+1) * n / 2
        ks1 = ks - n + 1
        gasdrag = - appu(ks1)
        ppx  = appu(ks) + appu(ks1)
        ppy  = appw(ks) + appw(ks1)
        ppdrag  = DSQRT( ppx**2 + ppy**2 ) 
        gas2pp(is,ijk) = ppdrag / gasdrag
        IF ( rlk(is) == 0.D0 ) gas2pp(is,ijk) = 0.D0
        IF ( (ppdrag < 0.D0) .OR. (gasdrag < 0.D0) ) THEN
          CALL error('drag_ratio','control gasdrag and ppdrag',1)
        END IF
      END DO
!
      RETURN
      END SUBROUTINE drag_ratio
!------------------------------------------------------------------------
      SUBROUTINE grav_ratio(kpgv, rlk, wpart, wgas, ijk)
! ... compute the ratio between gas-particle drag and the sum of particle-particle 
! ... drag coefficients for each solid phase
!
      USE atmosphere, ONLY: gravz
      USE time_parameters, ONLY: time, dt
      IMPLICIT NONE

      REAL*8,  INTENT(IN) :: kpgv(:)
      REAL*8,  INTENT(IN) :: rlk(:), wpart(:)
      REAL*8,  INTENT(IN) :: wgas
      INTEGER, INTENT(IN) :: ijk
      INTEGER :: is, n, ks, ks1
      REAL*8 :: gravterm, dragterm
!
      DO is = 1, nsolid
        gravterm = - gravz * rlk(is)
        dragterm = kpgv(is) * DABS(wpart(is) - wgas)
        drag2grav(is,ijk) = dragterm / gravterm
        IF ( rlk(is) == 0.D0 ) drag2grav(is,ijk) = 0.D0
      END DO
!
      RETURN
      END SUBROUTINE grav_ratio
!----------------------------------------------------------------------
      SUBROUTINE richardson(p,rgp,rhog,rlk,ug,wg,uk,wk,ijk,ijpk)
!
        USE atmosphere, ONLY : gravz, p_atm, t_atm
        USE domain_decomposition, ONLY: myijk
        USE gas_constants, ONLY: gmw, rgas, gammaair, gamn
        USE grid, ONLY: zb, dz
        USE indijk_module, ONLY: ip0_jp0_kp0_
        USE particles_constants, ONLY: inrl
        IMPLICIT NONE
!
        REAL*8,  INTENT(IN) :: p(:)
        REAL*8,  INTENT(IN) :: rgp(:)
        REAL*8,  INTENT(IN) :: rhog
        REAL*8,  INTENT(IN) :: rlk(:,:)
        REAL*8,  INTENT(IN) :: ug,wg
        REAL*8,  INTENT(IN) :: uk(:),wk(:)
        INTEGER, INTENT(IN) :: ijk, ijpk
!
        INTEGER :: j, imesh, is
        REAL*8 :: rhok, rhom, rhom1, um, wm
        REAL*8 :: p0, p1
        REAL*8 :: rhorif, prif, trif
        REAL*8 :: zrif
        REAL*8 :: epsk, ingam, factor
!
        imesh = myijk( ip0_jp0_kp0_, ijk)
        j = ( imesh - 1 ) / nx + 1

        zrif=zb(j)+0.5D0*(dz(1)-dz(j))
        prif = p_atm(j)
        trif = t_atm(j)

        rhorif = prif*gmw(6)/(rgas*trif)

        p0 = p(ijk)
        p1 = p(ijpk)
        rich(ijk) = 0.D0
        epsk = 0.D0
        rhom = 0.D0
        rhom1 = 0.D0
        um   = rgp(ijk) * ug
        wm   = rgp(ijk) * wg
        DO is = 1, nsolid
          epsk = epsk + rlk(ijk,is)*inrl(is)
          rhom = rhom + rlk(ijk,is)
          rhom1 = rhom1 + rlk(ijpk,is)
          um = um + rlk(ijk,is) * uk(is)
          wm = wm + rlk(ijk,is) * wk(is)
        END DO
        rhok = rhom
        rhom = rhom + rgp(ijk) 
        rhom1 = rhom1 + rgp(ijpk) 
        um = um / rhom
        wm = wm / rhom
 
        ingam = 1.D0 / gammaair
        factor = p0**ingam/(rhog+rhok)/gamn
        IF (epsk > 1E-8 ) THEN
          rich(ijk) = factor * (p1**gamn - p0**gamn) - gravz * 0.5*(dz(j)+dz(j+1))
          rich(ijk) = rich(ijk) / (um**2+wm**2)
          rich(ijk) = - ( rhom - rhom1 ) * gravz * dz(j) / ( rhom*(um**2+wm**2) )
        ELSE 
          rich(ijk) = 0.D0
        END IF
 
        RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE print_numbers
      USE time_parameters, ONLY : time
      IMPLICIT NONE

      INTEGER :: ij2, ij1, is, j
      INTEGER :: ijl
      INTEGER :: ig

      OPEN(UNIT=3,FILE='OUTPUT_NUMBERS')

      WRITE(3,547) time

      WRITE(3,101)
      DO j=1,nz
        ij1=1+(nz-j)*nx
        ij2=nx+(nz-j)*nx
        WRITE(3,550)(rich(ijl),ijl=ij1,ij2)
      END DO

      WRITE(3,102)
      DO j=1,nz
        ij1=1+(nz-j)*nx
        ij2=nx+(nz-j)*nx
        WRITE(3,550)(mut2mu(ijl),ijl=ij1,ij2)
      END DO

      DO is = 1, nsolid
        WRITE(3,103)is
        DO j=1,nz
          ij1=1+(nz-j)*nx
          ij2=nx+(nz-j)*nx
          WRITE(3,550)(gas2pp(is,ijl),ijl=ij1,ij2)
        END DO
      END DO

      DO is = 1, nsolid
        WRITE(3,104)is
        DO j=1,nz
          ij1=1+(nz-j)*nx
          ij2=nx+(nz-j)*nx
          WRITE(3,550)(drag2grav(is,ijl),ijl=ij1,ij2)
        END DO
      END DO

      CLOSE (3)
!
 547  FORMAT(1x,///,1x,'@@@ TIME = ',G11.4)
 550  FORMAT(1x,10(1x,g14.6e3))
 101  FORMAT(1x,//,1x,'RICHARDSON',/)
 102  FORMAT(1x,//,1x,'MUGT/MUG',/)
 103  FORMAT(1x,//,1x,'DRAG/PP',I1,/)
 104  FORMAT(1x,//,1x,'GRAV/DRAG',I1,/)
      END SUBROUTINE
!
!----------------------------------------------------------------------
      END MODULE nondim_numbers
!------------------------------------------------------------------------
