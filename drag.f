!----------------------------------------------------------------------
      MODULE momentum_transfer
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE kdrags(kpgv,du,dv,dw,ep,rgp,rlk,mug,is)
!----------------------------------------------------------------------
! ... This routine computes the gas-particles drag coefficient 
!
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: kpgv
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
      kpgv = drag
!
      RETURN
      END SUBROUTINE
!
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
      SUBROUTINE kdragg(drag,vrel, ep, rgp, rlk, mug, is)
!----------------------------------------------------------------------
      USE particles_constants, ONLY: rl, inrl, phis, dk
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: drag
      REAL*8, INTENT(IN) :: vrel, ep, rgp, rlk, mug
      INTEGER, INTENT(IN) :: is
      REAL*8 :: reynum, drvrel, dracoe
!
      dracoe=4.4D-1
      reynum = rgp * dk(is) * phis(is) * vrel/mug
! ... ?
      IF(reynum < 1.D-3) reynum=1.D-3
!
      IF(reynum <= 1000.D0) THEN
        dracoe = (24.D0/reynum) * (1.D0 + 0.15D0*reynum**0.687D0)
      END IF

      drvrel = dracoe * vrel / ep**2.7D0

      IF(drvrel > 1.D30) THEN
        drag=1.D30
      ELSE
        drag = 0.75D0 * drvrel * rgp / (dk(is)*phis(is))
        drag = drag * rlk * inrl(is) 
      END IF 
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE kdragl(drag, vrel, ep, rgp, rlk, mug, is)
!
      USE particles_constants, ONLY: rl, inrl, phis, dk
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: drag
      REAL*8, INTENT(IN) :: vrel, ep, rgp, rlk, mug
      INTEGER, INTENT(IN) :: is
      REAL*8 :: denom
!
      denom = 1.0D0 / (dk(is) * phis(is) * ep)
      drag = 150.D0 * ep * rlk * inrl(is) * mug  * denom +   &
             1.75D0 * rgp * vrel
      drag = drag  * rlk * inrl(is) * denom
      CONTINUE
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
        !effe=(3.D0*epkl**(1.D0/3.D0)+epsum**(1.D0/3.D0))/ &
        !     (2.D0*(epkl**(1.D0/3.D0)-epsum**(1.D0/3.D0)))
        effe = (3.D0*m1+m2) / (2.D0*(m1-m2))
        ppdr = fac*(1.D0+restc)*rlk(k-1)*rlk(kk-1)*dkf(k-1,kk-1)*effe
      ELSE
        ppdr = 0.D0
      ENDIF
      RETURN
      END SUBROUTINE
!------------------------------------------------------------------------
      END MODULE momentum_transfer
!------------------------------------------------------------------------
