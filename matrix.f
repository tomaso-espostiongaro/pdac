!--------------------------------------------------------------------
      MODULE phases_matrix
!--------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
      REAL*8, DIMENSION(:),   ALLOCATABLE ::  bu1, bv1, bu, bv
      REAL*8, DIMENSION(:,:), ALLOCATABLE ::  au1, av1, au, av
      REAL*8 :: rlim
!--------------------------------------------------------------------
      CONTAINS
!--------------------------------------------------------------------
      SUBROUTINE bounds_matrix
      USE dimensions
      IMPLICIT NONE
!
      ALLOCATE(bu1(nphase), bv1(nphase), bu(nphase), bv(nphase))
      ALLOCATE(au1(nphase,nphase), av1(nphase,nphase),           &
                au(nphase,nphase),  av(nphase,nphase))
!
      au = 0.D0;   av = 0.D0
      au1 = 0.D0;  av1 = 0.D0
      bu = 0.D0;   bv = 0.D0
      bu1 = 0.D0;  bv1 = 0.D0
!
      RETURN
      END SUBROUTINE

! ... MODIFICARE_X3D ( fino fine file )

!--------------------------------------------------------------------
      SUBROUTINE mats(ij, ep, epr, ept, epl, epb, p, pr, pt, pl, pb,  &
        rlk, rlkr, rlkt, rlkl, rlkb, rgp, rgpr, rgpt, rgpl, rgpb)
!
! ... Computes matrix elements to solve momentum-balance 
! ... linear system of coupled equations, in current cell and in
! ... left and bottom cells 
!
      USE dimensions
      USE grid, ONLY: dz, dr
      USE grid, ONLY: myij, myinds
      USE grid, ONLY: r_, t_, l_, b_
      USE tilde_momentum, ONLY: appu, appv
      USE particles_constants, ONLY: rl, inrl
      USE tilde_momentum, ONLY: rug, rvg, ruk, rvk
      USE time_parameters, ONLY: dt

      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ij
!
      INTEGER :: ks, kk, ks1, k
      INTEGER :: i,j,imesh
      INTEGER :: ijm, imj, ijr, ijt, ijl, ijb
      REAL*8 :: drp, dzp, dzm, drm, indrp, indzp, indzm, indrm
      REAL*8 :: ep_w, ep_s, ep_e, ep_n
      REAL*8 :: epk_w, epk_s, epk_e, epk_n
      REAL*8, INTENT(IN) :: ep, epr, ept, epl, epb, p, pr, pt, pl, pb
      REAL*8, INTENT(IN) :: rlk(:), rlkr(:), rlkt(:), rlkl(:), rlkb(:)
      REAL*8, INTENT(IN) :: rgp, rgpr, rgpt, rgpl, rgpb
!
      imesh = myij(0, 0, ij)
      j  = ( imesh - 1 ) / nr + 1
      i  = MOD( ( imesh - 1 ), nr) + 1
      ijm = myij( 0,-1, ij)
      imj = myij(-1, 0, ij)
      ijr = myinds(r_, ij)
      ijt = myinds(t_, ij)
      ijl = myinds(l_, ij)
      ijb = myinds(b_, ij)
!
      drm=dr(i)+dr(i-1)
      dzm=dz(j)+dz(j-1)
!
      indrm=1.D0/drm
      indzm=1.D0/dzm
!
! ... Backward ...
!
      DO k=1,nphase
!
! ... Explicit terms in the linear system
!
        IF(k.EQ.1) THEN
          ep_w = (dr(i)*epl + dr(i-1)*ep) * indrm
          ep_s = (dz(j)*epb + dz(j-1)*ep) * indzm
          bu1(1) = rug(imj)+ dt * indrm *2.D0* ep_w * (pl-p)
          bv1(1) = rvg(ijm)+ dt * indzm *2.D0* ep_s * (pb-p)
        ELSE
          epk_w = (dr(i)*rlkl(k-1) + dr(i-1)*rlk(k-1)) * indrm * inrl(k-1)
          epk_s = (dz(j)*rlkb(k-1) + dz(j-1)*rlk(k-1)) * indzm * inrl(k-1)
          bu1(k) = ruk(k-1,imj) + dt * indrm *2.D0* epk_w * (pl-p)
          bv1(k) = rvk(k-1,ijm) + dt * indzm *2.D0* epk_s * (pb-p)
        ENDIF
!
! ... Implicit terms in the linear system
!
        ks1=k*(k-1)/2
        DO kk=1,k
          ks=ks1+kk
          au1(k,kk)=(dr(i)*appu(ks,ijl)+appu(ks,ij)*dr(i-1))*indrm
          au1(kk,k)=au1(k,kk)
          av1(k,kk)=(dz(j)*appv(ks,ijb)+dz(j-1)*appv(ks,ij))*indzm
          av1(kk,k)=av1(k,kk)
        END DO
        IF(k.EQ.1) THEN
          au1(1,1)=au1(1,1)+(dr(i)*rgpl+dr(i-1)*rgp)*indrm
          av1(1,1)=av1(1,1)+(dz(j)*rgpb+dz(j-1)*rgp)*indzm
        ELSE
          au1(k,k)=au1(k,k)+(dr(i)*rlkl(k-1)+dr(i-1)*rlk(k-1))*indrm
          av1(k,k)=av1(k,k)+(rlkb(k-1)*dz(j)+dz(j-1)*rlk(k-1))*indzm
        ENDIF
      END DO
!
! ... Forward ...
!
      drp=dr(i)+dr(i+1)
      dzp=dz(j)+dz(j+1)
!
      indrp=1.D0/drp
      indzp=1.D0/dzp
!
      DO k=1,nphase
!
! ... Explicit terms in the linear system
!
        IF(k.EQ.1) THEN
          ep_e = (dr(i)*epr + dr(i+1)*ep) * indrp
          ep_n = (dz(j)*ept + dz(j+1)*ep) * indzp
          bu(1)  = rug(ij)+ dt * indrp *2.D0* ep_e * (p-pr)
          bv(1)  = rvg(ij)+ dt * indzp *2.D0* ep_n * (p-pt)
        ELSE
          epk_e = (dr(i)*rlkr(k-1) + dr(i+1)*rlk(k-1)) * indrp * inrl(k-1)
          epk_n = (dz(j)*rlkt(k-1) + dz(j+1)*rlk(k-1)) * indzp * inrl(k-1)
          bu(k)  = ruk(k-1,ij) + dt * indrp *2.D0* epk_e * (p-pr)
          bv(k)  = rvk(k-1,ij) + dt * indzp *2.D0* epk_n * (p-pt)
        ENDIF
!
! ... Implicit terms in the linear system
!
        ks1=k*(k-1)/2
        DO kk=1,k
          ks=ks1+kk
          au(k,kk)=(dr(i)*appu(ks,ijr)+appu(ks,ij)*dr(i+1))*indrp
          au(kk,k)=au(k,kk)
          av(k,kk)=(dz(j)*appv(ks,ijt)+dz(j+1)*appv(ks,ij))*indzp
          av(kk,k)=av(k,kk)
        END DO
        IF(k.EQ.1) THEN
          au(1,1)=au(1,1)+(dr(i)*rgpr+dr(i+1)*rgp)*indrp
          av(1,1)=av(1,1)+(dz(j)*rgpt+dz(j+1)*rgp)*indzp
        ELSE
          au(k,k)=au(k,k)+(dr(i)*rlkr(k-1)+dr(i+1)*rlk(k-1))*indrp
          av(k,k)=av(k,k)+(rlkt(k-1)*dz(j)+dz(j+1)*rlk(k-1))*indzp
        ENDIF
      END DO
!
      RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
      SUBROUTINE matsa(ij, ep, epr, ept, p, pr, pt, rlk, rlkr, rlkt, rgp, rgpr, rgpt)
!
      USE dimensions
      USE grid, ONLY: dz, dr
      USE grid, ONLY: myij, myinds
      USE grid, ONLY: r_, t_
      USE tilde_momentum, ONLY: appu, appv
      USE particles_constants, ONLY: rl, inrl
      USE tilde_momentum, ONLY: rug, rvg, ruk, rvk
      USE time_parameters, ONLY: dt

      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ij
      REAL*8, INTENT(IN) :: ep, epr, ept, p, pr, pt
      REAL*8, INTENT(IN) :: rlk(:), rlkr(:), rlkt(:)
      REAL*8, INTENT(IN) :: rgp, rgpr, rgpt
      INTEGER :: ks, kk, ks1, k
      INTEGER :: imesh, i, j
      INTEGER :: ijr, ijt
      REAL*8 :: drp, dzp, indrp, indzp
      REAL*8 :: ep_e, ep_n
      REAL*8 :: epk_e, epk_n
!
      imesh = myij( 0, 0, ij)
      j  = ( imesh - 1 ) / nr + 1
      i  = MOD( ( imesh - 1 ), nr) + 1
      ijr = myinds(r_, ij)
      ijt = myinds(t_, ij)
!
! ... Forward ...
!
      drp=dr(i)+dr(i+1)
      dzp=dz(j)+dz(j+1)
!
      indrp=1.D0/drp
      indzp=1.D0/dzp
!
      DO k=1,nphase
!
! ... Explicit terms in the linear system
!
        IF(k.EQ.1) THEN
          ep_e = (dr(i)*epr + dr(i+1)*ep) * indrp                  
          ep_n = (dz(j)*ept + dz(j+1)*ep) * indzp
          bu(1)  = rug(ij)+ dt * indrp *2.D0* ep_e * (p-pr)
          bv(1)  = rvg(ij)+ dt * indzp *2.D0* ep_n * (p-pt)
        ELSE
          epk_e = (dr(i)*rlkr(k-1) + dr(i+1)*rlk(k-1)) * indrp * inrl(k-1)
          epk_n = (dz(j)*rlkt(k-1) + dz(j+1)*rlk(k-1)) * indzp * inrl(k-1)
          bu(k)  = ruk(k-1,ij) + dt * indrp *2.D0* epk_e * (p-pr)
          bv(k)  = rvk(k-1,ij) + dt * indzp *2.D0* epk_n * (p-pt)
        ENDIF
!
! ... Implicit terms in the linear system
!
        ks1=k*(k-1)/2
        DO kk=1,k
          ks=ks1+kk
          au(k,kk)=(dr(i)*appu(ks,ijr)+appu(ks,ij)*dr(i+1))*indrp
          au(kk,k)=au(k,kk)
          av(k,kk)=(dz(j)*appv(ks,ijt)+dz(j+1)*appv(ks,ij))*indzp
          av(kk,k)=av(k,kk)
        END DO
        IF(k.EQ.1) THEN
          au(1,1)=au(1,1)+(dr(i)*rgpr+dr(i+1)*rgp)*indrp
          av(1,1)=av(1,1)+(dz(j)*rgpt+dz(j+1)*rgp)*indzp
        ELSE
          au(k,k)=au(k,k)+(dr(i)*rlkr(k-1)+dr(i+1)*rlk(k-1))*indrp
          av(k,k)=av(k,k)+(rlkt(k-1)*dz(j)+dz(j+1)*rlk(k-1))*indzp
        ENDIF
      END DO
!
      RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
      SUBROUTINE velsk(ug, vg, uk, vk, ugm, vgm, ukm, vkm, ij)
!
      USE grid, ONLY: fl_l
      USE grid, ONLY: myij
      USE dimensions
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: ug,vg,uk(:),vk(:),ugm,vgm,ukm(:),vkm(:)
      INTEGER, INTENT(IN) :: ij
!
      INTEGER :: kk, kp1, k, kj, ki
      INTEGER :: i,j,imesh
      INTEGER :: imj, ijm, ipj, ijp
      REAL*8 :: flt, flr, fll, div, amul, flb
!
      imesh = myij( 0, 0, ij)
      j  = ( imesh - 1 ) / nr + 1
      i  = MOD( ( imesh - 1 ), nr) + 1
!
! ... Use Gauss-Jordan method for matrix inversion
!
      imj = myij(-1, 0, ij)
      fll=fl_l(imj)
      IF(.NOT.(fll.EQ.2.OR.fll.eq.3.or.fll.eq.5)) THEN

        DO k=2,nphase
          IF(au1(k,k).LT.rlim) THEN
            DO kk=1,nphase
              au1(kk,kk)=au1(kk,kk)+au1(k,kk)
              au1(k,kk)=0.D0
              au1(kk,k)=0.D0
            END DO
            bu1(k)=0.D0
          END IF
        END DO

        DO k=1,nphase
          IF(au1(k,k).NE.0.D0) THEN 
            kp1=k+1
            div=1.D0/au1(k,k)
            DO kj=kp1,nphase
              au1(k,kj)=au1(k,kj)*div
            END DO
            bu1(k)=bu1(k)*div
            au1(k,k)=0.D0
            DO ki=1,nphase
              amul=au1(ki,k)
              DO kj=kp1,nphase
                au1(ki,kj)=au1(ki,kj)-amul*au1(k,kj)
              END DO
              bu1(ki)=bu1(ki)-amul*bu1(k)
            END DO
          END IF
        END DO
!
        ugm=bu1(1)
        DO k=2,nphase
          ukm(k-1)=bu1(k)
        END DO

      END IF

      ijm = myij( 0,-1, ij)
      flb=fl_l(ijm)
      IF(.NOT.(flb.EQ.2.OR.flb.eq.3.or.flb.eq.5)) THEN

        DO k=2,nphase
          IF(av1(k,k).LT.rlim) THEN
            DO kk=1,nphase
              av1(kk,kk)=av1(kk,kk)+av1(k,kk)
              av1(k,kk)=0.D0
              av1(kk,k)=0.D0
            END DO
            bv1(k)=0.D0
          END IF
        END DO

        DO k=1,nphase
          IF(av1(k,k).NE.0.D0) THEN
            kp1=k+1
            div=1.D0/av1(k,k)
            DO kj=kp1,nphase
              av1(k,kj)=av1(k,kj)*div
            END DO
            bv1(k)=bv1(k)*div
            av1(k,k)=0.D0
            DO ki=1,nphase
              amul=av1(ki,k)
              DO kj=kp1,nphase
                av1(ki,kj)=av1(ki,kj)-amul*av1(k,kj)
              END DO
              bv1(ki)=bv1(ki)-amul*bv1(k)
            END DO
          END IF
        END DO
!
        vgm=bv1(1)
        DO k=2, nphase
          vkm(k-1)=bv1(k)
        END DO
!
      END IF

      ipj = myij(+1, 0, ij)
      flr=fl_l(ipj)
      IF(.NOT.(flr.EQ.2.OR.flr.eq.3.or.flr.eq.5)) THEN

        DO k=2,nphase
          IF(au(k,k).LT.rlim) THEN
            DO kk=1,nphase
              au(kk,kk)=au(kk,kk)+au(k,kk)
              au(k,kk)=0.D0
              au(kk,k)=0.D0
            END DO
            bu(k)=0.D0
          END IF
        END DO

        DO k=1,nphase
          IF(au(k,k).NE.0.D0) THEN
            kp1=k+1
            div=1.D0/au(k,k)
            DO kj=kp1,nphase
              au(k,kj)=au(k,kj)*div
            END DO
            bu(k)=bu(k)*div
            au(k,k)=0.D0
            DO ki=1,nphase
              amul=au(ki,k)
              DO kj=kp1,nphase
                au(ki,kj)=au(ki,kj)-amul*au(k,kj)
              END DO
              bu(ki)=bu(ki)-amul*bu(k)
            END DO
          END IF
        END DO
!
        ug=bu(1)
        DO k=2,nphase
          uk(k-1)=bu(k)
        END DO
!
      END IF

      ijp = myij( 0,+1, ij)
      flt=fl_l(ijp)
      IF(.NOT.(flt.EQ.2.OR.flt.eq.3.or.flt.eq.5)) THEN

        DO k=2,nphase
          IF(av(k,k).LT.rlim) THEN
            DO kk=1,nphase
              av(kk,kk)=av(kk,kk)+av(k,kk)
              av(k,kk)=0.D0
              av(kk,k)=0.D0
            END DO
            bv(k)=0.D0
          END IF
        END DO

        DO k=1,nphase
          IF(av(k,k).NE.0.D0) THEN 
            kp1=k+1
            div=1.D0/av(k,k)
            DO kj=kp1,nphase
              av(k,kj)=av(k,kj)*div
            END DO
            bv(k)=bv(k)*div
            av(k,k)=0.D0
            DO ki=1,nphase
              amul=av(ki,k)
              DO kj=kp1,nphase
                av(ki,kj)=av(ki,kj)-amul*av(k,kj)
              END DO 
              bv(ki)=bv(ki)-amul*bv(k)
            END DO 
          END IF
        END DO
!
        vg=bv(1)
        DO k=2,nphase
          vk(k-1)=bv(k)
        END DO
!
      END IF
!
      RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
      SUBROUTINE velsk2(ug, vg, uk, vk, ij)
!
      USE grid, ONLY: fl_l
      USE grid, ONLY: myij
      USE dimensions
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: ug, vg, uk(:), vk(:)
      INTEGER, INTENT(IN) :: ij
      INTEGER :: ipj, ijp
!
      INTEGER :: kk, kp1, k, kj, ki
      REAL*8 :: flt, flr, fll, div, amul, flb

      ipj = myij(+1, 0, ij)
      flr=fl_l(ipj)
      IF(.NOT.(flr.EQ.2.OR.flr.eq.3.or.flr.eq.5)) THEN

        DO k=2,nphase
          IF(au(k,k).LT.rlim) THEN
            DO kk=1,nphase
              au(kk,kk)=au(kk,kk)+au(k,kk)
              au(k,kk)=0.D0
              au(kk,k)=0.D0
            END DO
            bu(k)=0.D0
          END IF
        END DO

        DO k=1,nphase
          IF(au(k,k).NE.0.D0) THEN
            kp1=k+1
            div=1.D0/au(k,k)
            DO kj=kp1,nphase
              au(k,kj)=au(k,kj)*div
            END DO
            bu(k)=bu(k)*div
            au(k,k)=0.D0
            DO ki=1,nphase
              amul=au(ki,k)
              DO kj=kp1,nphase
                au(ki,kj)=au(ki,kj)-amul*au(k,kj)
              END DO
              bu(ki)=bu(ki)-amul*bu(k)
            END DO
          END IF
        END DO
!
        ug=bu1(1)
        DO k=2,nphase
          uk(k-1)=bu(k)
        END DO
!
      END IF

      ijp = myij( 0,+1, ij)
      flt=fl_l(ijp)
      IF(.NOT.(flt.EQ.2.OR.flt.eq.3.or.flt.eq.5)) THEN

        DO k=2,nphase
          IF(av(k,k).LT.rlim) THEN
            DO kk=1,nphase
              av(kk,kk)=av(kk,kk)+av(k,kk)
              av(k,kk)=0.D0
              av(kk,k)=0.D0
            END DO
            bv(k)=0.D0
          END IF
        END DO

        DO k=1,nphase
          IF(av(k,k).NE.0.D0) THEN 
            kp1=k+1
            div=1.D0/av(k,k)
            DO kj=kp1,nphase
              av(k,kj)=av(k,kj)*div
            END DO
            bv(k)=bv(k)*div
            av(k,k)=0.D0
            DO ki=1,nphase
              amul=av(ki,k)
              DO kj=kp1,nphase
                av(ki,kj)=av(ki,kj)-amul*av(k,kj)
              END DO 
              bv(ki)=bv(ki)-amul*bv(k)
            END DO 
          END IF
        END DO
!
        vg=bv(1)
        DO k=2,nphase
          vk(k-1)=bv(k)
        END DO
!
      END IF
!
      RETURN
      END SUBROUTINE

!------------------------------------------------------------------------
      END MODULE
!------------------------------------------------------------------------
