!------------------------------------------------------------------------
      MODULE nondim_numbers
!------------------------------------------------------------------------
      USE dimensions
      USE grid, ONLY: ncint
      IMPLICIT NONE
      SAVE
!
      REAL*8, ALLOCATABLE :: ratio(:,:)   ! Gas-Particle / Particle-Particle ratio !
      REAL*8, ALLOCATABLE :: rich(:)      ! Richardson number
      REAL*8, ALLOCATABLE :: mut2mu(:)    ! Turbulent/Molecular viscosity
!------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------
      SUBROUTINE drag_ratio(appu,appw,kpgv,ijk)
! ... compute the ratio between gas-particle drag and the sum of particle-particle 
! ... drag coefficients for each solid phase
!
      USE time_parameters, ONLY: time, dt
      USE gas_solid_density, ONLY : rlk
      IMPLICIT NONE

      REAL*8,  INTENT(IN) :: appu(:), appw(:)
      REAL*8,  INTENT(IN) :: kpgv(:)
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
        ratio(is,ijk) = ppdrag / gasdrag
        IF ( rlk(is,ijk) == 0.D0 ) ratio(is,ijk) = 0.D0
        IF ( (ppdrag < 0.D0) .OR. (gasdrag < 0.D0) ) THEN
          CALL error('drag_ratio','control gasdrag and ppdrag',1)
        END IF
      END DO
!
      RETURN
      END SUBROUTINE drag_ratio
!----------------------------------------------------------------------
      SUBROUTINE richardson(rgp,rlk,ug,wg,uk,wk,ijk,ijpk)
!
        USE atmosphere, ONLY : gravz, atm
        USE gas_constants, ONLY: gmw, rgas
        USE particles_constants, ONLY: inrl
        USE grid, ONLY: myijk, zb, dz
        USE indijk_module, ONLY: ip0_jp0_kp0_
        IMPLICIT NONE
!
        REAL*8,  INTENT(IN) :: rgp(:)
        REAL*8,  INTENT(IN) :: rlk(:,:)
        REAL*8,  INTENT(IN) :: ug,wg
        REAL*8,  INTENT(IN) :: uk(:),wk(:)
        INTEGER, INTENT(IN) :: ijk, ijpk
!
        INTEGER :: j, imesh, is
        REAL*8 :: rhom, rhom1, um, wm
        REAL*8 :: rhorif, prif, trif
        REAL*8 :: zrif
        REAL*8 :: epsk
!
        imesh = myijk( ip0_jp0_kp0_, ijk)
        j = ( imesh - 1 ) / nr + 1
        zrif=zb(j)+0.5D0*(dz(1)-dz(j))

        CALL atm(zrif,prif,trif)
        rhorif = prif*gmw(6)/(rgas*trif)

        epsk = 0.D0
        rich(ijk) = 0.D0
        rhom = rgp(ijk) 
        rhom1 = rgp(ijpk) 
        um   = rgp(ijk) * ug
        wm   = rgp(ijk) * wg
        DO is = 1, nsolid
          epsk = epsk + rlk(is,ijk)*inrl(is)
          rhom = rhom + rlk(is,ijk)
          rhom1 = rhom1 + rlk(is,ijpk)
          um = um + rlk(is,ijk) * uk(is)
          wm = wm + rlk(is,ijk) * wk(is)
        END DO
        um = um / rhom
        wm = wm / rhom

        IF (epsk > 1E-6 ) THEN
          rich(ijk) = ((rhom - rhom1) * gravz * dz(j))/(rhom*(um**2+wm**2))
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
        ij1=1+(nz-j)*nr
        ij2=nr+(nz-j)*nr
        WRITE(3,550)(rich(ijl),ijl=ij1,ij2)
      END DO

      WRITE(3,103)
      DO j=1,nz
        ij1=1+(nz-j)*nr
        ij2=nr+(nz-j)*nr
        WRITE(3,550)(mut2mu(ijl),ijl=ij1,ij2)
      END DO

      DO is = 1, nsolid
        WRITE(3,102)is
        DO j=1,nz
          ij1=1+(nz-j)*nr
          ij2=nr+(nz-j)*nr
          WRITE(3,550)(ratio(is,ijl),ijl=ij1,ij2)
        END DO
      END DO

      CLOSE (3)
!
 547  FORMAT(1x,///,1x,'@@@ TIME = ',G11.4)
 550  FORMAT(1x,10(1x,g14.6e3))
 101  FORMAT(1x,//,1x,'RICHARDSON',/)
 103  FORMAT(1x,//,1x,'MUGT/MUG',/)
 102  FORMAT(1x,//,1x,'DRAG/PP',I1,/)
      END SUBROUTINE
!
!----------------------------------------------------------------------
      END MODULE nondim_numbers
!------------------------------------------------------------------------
