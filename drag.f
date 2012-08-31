!----------------------------------------------------------------------
      MODULE momentum_transfer
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: drag_model
      SAVE
      PRIVATE :: kdragl, kdragg, f, ppdrag
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE sdrag(drag,du,dv,dw,ep,rgp,rlk,mug,is)
!----------------------------------------------------------------------
! ... This routine computes the gas-particles drag coefficient using Syamlal and O'Brien
! ... continuous correlation
!
      USE particles_constants, ONLY: rl, inrl, phis, dk
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) ::  du, dv, dw
      REAL*8, INTENT(IN) :: ep,rgp,mug
      REAL*8, INTENT(IN) :: rlk
      INTEGER, INTENT(IN) :: is
      REAL*8 :: drag, vrel, drag_coeff
      REAL*8 :: reynum, eps
      REAL*8 :: a, b, vr
!
      vrel = DSQRT(du**2 + dv**2 + dw**2)
      eps = rlk * inrl(is)
      reynum = rgp * dk(is) * phis(is) * vrel/mug
!
! ... Ratio of the terminal settling velocity of a multiparticle system to that
! ... of an isolated single particle
!
      a = ep**4.14
      IF (ep <= 0.85) THEN
        b = 0.8*ep**1.28
      ELSE
        b = ep**2.65
      END IF
!
      vr = 0.5D0 * (a - 0.06*reynum + &
           DSQRT(0.0036*reynum**2 + 0.12*reynum*(2*b-a) + a**2))
!
! ... Drag coefficient (Dalla Valle, 1948)
!
      drag_coeff = (0.63D0 + 4.8D0/DSQRT(reynum))**2
!
! ... This is needed for convergence at low reynum
      drag_coeff = MIN(drag_coeff, 1.D4) 
!
! ... Fluid-particle drag (single particle)
!
      drag = 0.75D0 * drag_coeff * eps / (dk(is)*phis(is)) * rgp * vrel 
!
      drag = drag / vr
!
      RETURN
      END SUBROUTINE sdrag
!----------------------------------------------------------------------
      SUBROUTINE kdrags(drag,du,dv,dw,ep,rgp,rlk,mug,is)
!----------------------------------------------------------------------
! ... This routine computes the gas-particles drag coefficient 
!
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) ::  du, dv, dw
      REAL*8, INTENT(IN) :: ep,rgp,mug
      REAL*8, INTENT(IN) :: rlk
      INTEGER, INTENT(IN) :: is
      REAL*8 :: drag, vrel
!
      vrel = DSQRT(du**2 + dv**2 + dw**2)

      IF(ep < 0.8D0) THEN
        CALL kdragl(drag,vrel,ep,rgp,rlk,mug,is)
      ELSE
        CALL kdragg(drag,vrel,ep,rgp,rlk,mug,is)
      END IF
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE kdragg(drag,vrel, ep, rgp, rlk, mug, is)
!
! ... Wen and Yu (1966) expression for dilute regime (ep > 0.8)
!
      USE particles_constants, ONLY: rl, inrl, phis, dk, phi
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: drag
      REAL*8, INTENT(IN) :: vrel, ep, rgp, rlk, mug
      INTEGER, INTENT(IN) :: is
      REAL*8 :: reynum, drag_coeff, eps
!
      eps = rlk * inrl(is)
      reynum = rgp * dk(is) * phis(is) * vrel/mug
!
! ... Drag coefficient
!
      IF(reynum <= 1000.D0) THEN
        ! ... modified Stokes law
        drag_coeff = (24.D0/reynum) * (1.D0 + 0.15*reynum**0.687)
      ELSE
        drag_coeff = 0.43828814
      END IF
      drag_coeff = MIN(drag_coeff, 1.D4) 
!
! ... Fluid-particle drag (single particle)
!
      drag = 0.75D0 * drag_coeff * eps / (dk(is)*phis(is)) * rgp * vrel 
!
! ... correction for mixtures
!
      drag = drag * f(ep)
!      drag = drag * fs(eps,phi(is))
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE kdragl(drag, vrel, ep, rgp, rlk, mug, is)
!
! ... The famous Ergun [1952] equation for dense regimes
!
      USE particles_constants, ONLY: rl, inrl, phis, dk
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: drag
      REAL*8, INTENT(IN) :: vrel, ep, rgp, rlk, mug
      INTEGER, INTENT(IN) :: is
      REAL*8 :: denom
!
! ... Gidaspow [1994] pag. 36
!
      denom = 1.0D0 / (dk(is) * phis(is) * ep)
!
      drag = 150.D0 * ep * rlk * inrl(is) * mug  * denom +   &
             1.75D0 * rgp * vrel
!
      drag = drag  * rlk * inrl(is) * denom
      CONTINUE
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      REAL*8 FUNCTION f(void_fraction)
!
! ... Ratio of fall velocity of a single particle with respect to a cloud
! ... Modification of the Richardson and Zaki (1954) equation
! ... [Gidaspow, 1994; Chapt.2]
!
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: void_fraction
      !
      f = void_fraction**(-2.65)
      !
      RETURN
      END FUNCTION f
!----------------------------------------------------------------------
      REAL*8 FUNCTION fs(solid_fraction,phimax)
!
! ... Ratio of fall velocity of a single particle with respect to a cloud
! ... Modification of the Richardson and Zaki (1954) equation
! ... [Michaels and Bolger, 1962]
!
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: solid_fraction, phimax
      !
      fs = (1.D0 - solid_fraction/phimax)**(-4.65)
      !
      RETURN
      END FUNCTION fs
!----------------------------------------------------------------------
      SUBROUTINE inter(appu,appv,appw,kpgv,us,vs,ws,rlk,ijk)
!----------------------------------------------------------------------
! ... This routine computes the interphase terms
!
      USE dimensions
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE time_parameters, ONLY: dt
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: appu(:), appv(:), appw(:)
      REAL*8, INTENT(IN) :: kpgv(:) 
      REAL*8, INTENT(IN) :: us(:,:),vs(:,:),ws(:,:)
      REAL*8, INTENT(IN) :: rlk(:,:)
      INTEGER, INTENT(IN) :: ijk
!
      INTEGER :: k, ks, kk
      REAL*8  :: sumx, sumy, sumz
      REAL*8 :: con
      REAL*8 :: du, dv, dw
!
      ks=0
      DO k=1,nphase
        IF(k == 1) THEN
          ks=ks+1
          appu(ks)=0.D0
          appv(ks)=0.D0
          appw(ks)=0.D0
        ELSE
          DO kk=1,k
              ks=ks+1
            IF(kk == 1)THEN
              appu(ks)= - kpgv(k-1) * dt
              appv(ks)= - kpgv(k-1) * dt
              appw(ks)= - kpgv(k-1) * dt
            ELSEIF(kk == k)THEN
              appu(ks)=0.D0
              appv(ks)=0.D0
              appw(ks)=0.D0
            ELSE
              CALL ppdrag(con,rlk(ijk,:),k,kk)
              du=DABS(us(ijk,k-1)-us(ijk,kk-1)+us(imjk,k-1)-us(imjk,kk-1))*0.5D0
              dv=DABS(vs(ijk,k-1)-vs(ijk,kk-1)+vs(ijmk,k-1)-vs(ijmk,kk-1))*0.5D0
              dw=DABS(ws(ijk,k-1)-ws(ijk,kk-1)+ws(ijkm,k-1)-ws(ijkm,kk-1))*0.5D0
              appu(ks)= - (con*du) * dt
              appv(ks)= - (con*dv) * dt
              appw(ks)= - (con*dw) * dt
            ENDIF
          END DO
        ENDIF
      END DO
      DO k=1,nphase
        sumx=0.D0
        sumy=0.D0
        sumz=0.D0
        DO kk=1,k
          ks=kk+k*(k-1)/2
          sumx=sumx+appu(ks)
          sumy=sumy+appv(ks)
          sumz=sumz+appw(ks)
        END DO
        DO kk=k,nphase
          ks=k+kk*(kk-1)/2
          sumx=sumx+appu(ks)
          sumy=sumy+appv(ks)
          sumz=sumz+appw(ks)
        END DO
        ks=k*(k+1)/2
        appu(ks)=-sumx
        appv(ks)=-sumy
        appw(ks)=-sumz
      END DO
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE ppdrag(ppdr,rlk,k,kk)
!
! ... This routine computes particle-particle drag coefficient
!
      USE particles_constants, ONLY: phi, epsl, dkf, epsu, dk, rl, inrl, philim
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: ppdr
      REAL*8, INTENT(IN) :: rlk(:)
      INTEGER, INTENT(IN) :: k, kk
!
      REAL*8 :: ep1, ep2, epsum
      REAL*8 :: xbar, epkl, effe
      REAL*8 :: fac, restc
      REAL*8 :: espo, m1, m2 
      INTEGER :: k1, k2
!
      fac=1.D0
      restc=1.D0
!
      IF(dk(k-1) >= dk(kk-1)) THEN
        k1=k-1
        k2=kk-1
      ELSE
        k1=kk-1
        k2=k-1
      ENDIF
      ep1=rlk(k1)*inrl(k1)
      ep2=rlk(k2)*inrl(k2)
      epsum=ep1+ep2
      IF(ep1 > 0.D0 .AND. ep2 > 0.D0) THEN
        xbar=ep1/epsum
        IF(xbar <= philim(k-1,kk-1)) THEN
          epkl = epsl(k-1,kk-1)*xbar + phi(k2)
        ELSE
          epkl = epsu(k-1,kk-1)*(1.D0-xbar) + phi(k1)
        ENDIF
!
! ... for Nakamura & Capes correlation effe=1.5
        !effe=1.5D0
        espo = 1.D0/3.D0  
        m1 = epkl**espo
        m2 = epsum**espo
        effe = ( 3.D0*m1 + m2 ) / ( 2.D0*( m1 - m2 ) )
        ppdr = fac*(1.D0+restc)*rlk(k-1)*rlk(kk-1)*dkf(k-1,kk-1)*effe
      ELSE
        ppdr = 0.D0
      ENDIF
              
      RETURN
      END SUBROUTINE ppdrag
!------------------------------------------------------------------------
      END MODULE momentum_transfer
!------------------------------------------------------------------------
