!----------------------------------------------------------------------
      MODULE momentum_transfer
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE kdrags(kpgv,ug,ugm,vg,vgm,uk,ukm,vk,vkm,ep,rog,rgp,rlk,mug)
!----------------------------------------------------------------------
! ... This routine computes the gas-particles drag coefficient 
!
      USE dimensions
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: kpgv(:)
      REAL*8, INTENT(IN) :: ug,ugm,vg,vgm,uk(:),ukm(:),vk(:),vkm(:),ep
      REAL*8, INTENT(IN) :: rog,rgp,rlk(:)
      REAL*8, INTENT(IN) :: mug
      INTEGER :: k
      REAL*8 :: dv, drag, dvgs, dugs
!
!..............................................................
! ... (legend)
!
!     X == X(ij)
!     ugm == ug(imj)
!     vgm == vg(ijm)
!     ukm(k) == uk(k,imj)
!     vkm(k) == vk(k,ijm)
!..............................................................
!
      DO k=1,ncl
        dvgs=DABS(( (vg-vk(k)) + (vgm-vkm(k)) ) * 0.5D0)
        dugs=DABS(( (ug-uk(k)) + (ugm-ukm(k)) ) * 0.5D0)
        dv=DSQRT(dvgs**2+dugs**2)
        IF(ep.LT.0.8D0) THEN
          CALL kdragl(drag,dv,ep,rog, rlk, k, mug)
        ELSE
          CALL kdragg(drag,dv,ep,rog, rgp, rlk, k, mug)
        END IF
        kpgv(k)=drag
      END DO
!
      RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
      SUBROUTINE inter(appu,appv,kpgv,uk,ukm,vk,vkm,rlk)
!----------------------------------------------------------------------
! ... This routine computes particle-particle drag coefficient
!
      USE dimensions
      USE particles_constants, ONLY: phi, epsl, dkf, epsu, dk, rl, inrl, philim
      USE time_parameters, ONLY: dt
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: appu(:), appv(:)
      REAL*8, INTENT(IN) :: kpgv(:) 
      REAL*8, INTENT(IN) :: uk(:),ukm(:),vk(:),vkm(:),rlk(:)
!
      INTEGER :: k, ks, kk, k1, k2
      REAL*8 :: xbar, epkl, ep2, epsum, sumx, sumy, effe, con
      REAL*8 :: ep1, fac, restc, dv, du
!
!..............................................................
! ... (legend)
!
!     X == X(ij)
!     ukm(k) == uk(k,imj)
!     vkm(k) == vk(k,ijm)
!..............................................................
!
      fac=1.D0
      restc=1.D0
!
      ks=0
      DO k=1,nphase
        IF(k.EQ.1) THEN
          ks=ks+1
          appu(ks)=0.D0
          appv(ks)=0.D0
        ELSE
          DO kk=1,k
              ks=ks+1
            IF(kk.EQ.1)THEN
              appu(ks)=-kpgv(k-1)*dt
              appv(ks)=-kpgv(k-1)*dt
            ELSEIF(kk.EQ.k)THEN
              appu(ks)=0.D0
              appv(ks)=0.D0
            ELSE
              dv=DABS(vk(k-1)-vk(kk-1) + vkm(k-1)-vkm(kk-1))*0.5D0
              du=DABS(uk(k-1)-uk(kk-1) + ukm(k-1)-ukm(kk-1))*0.5D0
              IF(dk(k-1).GE.dk(kk-1)) THEN
                k1=k-1
                k2=kk-1
              ELSE
                k1=kk-1
                k2=k-1
              ENDIF
              ep1=rlk(k1)*inrl(k1)
              ep2=rlk(k2)*inrl(k2)
              epsum=ep1+ep2
              IF(ep1.GT.0.D0.AND.ep2.gt.0.D0) THEN
                xbar=ep1/epsum
                IF(xbar.LE.philim(k-1,kk-1)) THEN
                  epkl=epsl(k-1,kk-1)*xbar/phi(k1)+phi(k2)
                ELSE
                  epkl=epsu(k-1,kk-1)*(1.D0-xbar)+phi(k1)
                ENDIF
!
! for Nakamura & Capes correlation effe=1.5
!                effe=1.5D0
                effe=(3.D0*epkl**(1.D0/3.D0)+epsum**(1.D0/3.D0))/    &
                     (2.D0*(epkl**(1.D0/3.D0)-epsum**(1.D0/3.D0)))
!
                con=fac*(1.D0+restc)*rlk(k-1)*rlk(kk-1)*dkf(k-1,kk-1)*dt*effe
              ELSE
                con=0.D0
              ENDIF
              appu(ks)=-con*du
              appv(ks)=-con*dv
            ENDIF
          END DO
        ENDIF
      END DO
      DO k=1,nphase
        sumx=0.D0
        sumy=0.D0
        DO kk=1,k
          ks=kk+k*(k-1)/2
          sumx=sumx+appu(ks)
          sumy=sumy+appv(ks)
        END DO
        DO kk=k,nphase
          ks=k+kk*(kk-1)/2
          sumx=sumx+appu(ks)
          sumy=sumy+appv(ks)
        END DO
        ks=k*(k+1)/2
        appu(ks)=-sumx
        appv(ks)=-sumy
      END DO
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE kdragg(drag,dv, ep, rog, rgp, rlk, k, mug)
!----------------------------------------------------------------------
      USE dimensions
      USE particles_constants, ONLY: rl, inrl, phis, dk
      IMPLICIT NONE
!
      INTEGER :: k
      REAL*8 :: reynum, drdvep, dracoe
      REAL*8, INTENT(OUT) :: drag
      REAL*8, INTENT(IN) :: dv, ep, rog, rgp, rlk(:), mug
!
      dracoe=4.4D-1
      reynum = rgp * dk(k) * phis(k) * dv/mug
      IF(reynum.LT.1.D-3) reynum=1.D-3

      IF(reynum.LE.1000.D0) THEN
        dracoe = (24.D0/reynum) * (1.D0 + 0.15D0*reynum**0.687D0)
      END IF

      drdvep = dracoe * dv / ep**2.7D0

      IF(drdvep.GT.1.D30) THEN
        drag=1.D30
      ELSE
        drag = 0.75D0 * drdvep * rog / (dk(k)*phis(k))
        drag = drag * rlk(k) * ep * inrl(k)
      END IF 
!
      RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
      SUBROUTINE kdragl(drag, dv, ep, rog, rlk, k, mug)
!
      USE dimensions
      USE particles_constants, ONLY: rl, inrl, phis, dk
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: k
      REAL*8 :: denom2, denom1
      REAL*8, INTENT(OUT) :: drag
      REAL*8, INTENT(IN) :: dv, ep, rog, rlk(:), mug
!
      denom2 = ep * dk(k) * phis(k)
      denom1 = denom2*denom2
      drag = 150.D0 * rlk(k) * inrl(k) * mug/denom1 + 1.75D0 * rog * dv/denom2
      drag = drag * rlk(k) * ep * inrl(k)
      CONTINUE
!
      RETURN
      END SUBROUTINE
!------------------------------------------------------------------------
      END MODULE momentum_transfer
!------------------------------------------------------------------------
