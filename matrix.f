!--------------------------------------------------------------------
      MODULE phases_matrix
!--------------------------------------------------------------------

      USE indijk_module

      IMPLICIT NONE
      SAVE

      REAL*8, PRIVATE, DIMENSION(:),   ALLOCATABLE ::  bu1, bw1, bu, bw
      REAL*8, PRIVATE, DIMENSION(:,:), ALLOCATABLE ::  au1, aw1, au, aw
      REAL*8 :: rlim
!--------------------------------------------------------------------
      CONTAINS
!--------------------------------------------------------------------
      SUBROUTINE bounds_matrix
      USE dimensions
      IMPLICIT NONE
!
      ALLOCATE(bu1(nphase), bw1(nphase), bu(nphase), bw(nphase))
      ALLOCATE(au1(nphase,nphase), aw1(nphase,nphase),           &
                au(nphase,nphase),  aw(nphase,nphase))
!
      au = 0.D0;   aw = 0.D0
      au1 = 0.D0;  aw1 = 0.D0
      bu = 0.D0;   bw = 0.D0
      bu1 = 0.D0;  bw1 = 0.D0
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
      USE grid, ONLY: myijk, myinds
      USE grid, ONLY: r_, t_, l_, b_
      USE tilde_momentum, ONLY: appu, appw
      USE particles_constants, ONLY: rl, inrl
      USE tilde_momentum, ONLY: rug, rwg, rus, rws
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
      REAL*8 :: eps_w, eps_s, eps_e, eps_n
      REAL*8, INTENT(IN) :: ep, epr, ept, epl, epb, p, pr, pt, pl, pb
      REAL*8, INTENT(IN) :: rlk(:), rlkr(:), rlkt(:), rlkl(:), rlkb(:)
      REAL*8, INTENT(IN) :: rgp, rgpr, rgpt, rgpl, rgpb
!
      imesh = myijk( ip0_jp0_kp0_, ij)
      j  = ( imesh - 1 ) / nr + 1
      i  = MOD( ( imesh - 1 ), nr) + 1
      ijm = myijk( ip0_jm1_kp0_, ij)
      imj = myijk( im1_jp0_kp0_, ij)
      ijr = myinds(ip1_jp0_kp0_, ij)
      ijt = myinds(ip0_jp1_kp0_, ij)
      ijl = myinds(im1_jp0_kp0_, ij)
      ijb = myinds(ip0_jm1_kp0_, ij)
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
          bw1(1) = rwg(ijm)+ dt * indzm *2.D0* ep_s * (pb-p)
        ELSE
          eps_w = (dr(i)*rlkl(k-1) + dr(i-1)*rlk(k-1)) * indrm * inrl(k-1)
          eps_s = (dz(j)*rlkb(k-1) + dz(j-1)*rlk(k-1)) * indzm * inrl(k-1)
          bu1(k) = rus(k-1,imj) + dt * indrm *2.D0* eps_w * (pl-p)
          bw1(k) = rws(k-1,ijm) + dt * indzm *2.D0* eps_s * (pb-p)
        ENDIF
!
! ... Implicit terms in the linear system
!
! ... number of elements of appu, appv, already used in previous cycles
! ... (= size of the upper triangular matrix of side k)
        ks1=k*(k-1)/2
!
! ... Off-diagonal elements
        DO kk=1,k
          ks=ks1+kk
          au1(k,kk)=(dr(i)*appu(ks,ijl)+appu(ks,ij)*dr(i-1))*indrm
          au1(kk,k)=au1(k,kk)
          aw1(k,kk)=(dz(j)*appw(ks,ijb)+dz(j-1)*appw(ks,ij))*indzm
          aw1(kk,k)=aw1(k,kk)
        END DO
!
! ... Diagonal elements
        IF(k.EQ.1) THEN
          au1(1,1)=au1(1,1)+(dr(i)*rgpl+dr(i-1)*rgp)*indrm
          aw1(1,1)=aw1(1,1)+(dz(j)*rgpb+dz(j-1)*rgp)*indzm
        ELSE
          au1(k,k)=au1(k,k)+(dr(i)*rlkl(k-1)+dr(i-1)*rlk(k-1))*indrm
          aw1(k,k)=aw1(k,k)+(rlkb(k-1)*dz(j)+dz(j-1)*rlk(k-1))*indzm
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
          bw(1)  = rwg(ij)+ dt * indzp *2.D0* ep_n * (p-pt)
        ELSE
          eps_e = (dr(i)*rlkr(k-1) + dr(i+1)*rlk(k-1)) * indrp * inrl(k-1)
          eps_n = (dz(j)*rlkt(k-1) + dz(j+1)*rlk(k-1)) * indzp * inrl(k-1)
          bu(k)  = rus(k-1,ij) + dt * indrp *2.D0* eps_e * (p-pr)
          bw(k)  = rws(k-1,ij) + dt * indzp *2.D0* eps_n * (p-pt)
        ENDIF
!
! ... Implicit terms in the linear system
!
! ... number of elements of appu, appv, already used in previous cycles
! ... (= size of the upper triangular matrix of side k)
        ks1=k*(k-1)/2
!
! ... Off-Diagonal elements
        DO kk=1,k
          ks=ks1+kk
          au(k,kk)=(dr(i)*appu(ks,ijr)+appu(ks,ij)*dr(i+1))*indrp
          au(kk,k)=au(k,kk)
          aw(k,kk)=(dz(j)*appw(ks,ijt)+dz(j+1)*appw(ks,ij))*indzp
          aw(kk,k)=aw(k,kk)
        END DO
!
! ... Diagonal elements
        IF(k.EQ.1) THEN
          au(1,1)=au(1,1)+(dr(i)*rgpr+dr(i+1)*rgp)*indrp
          aw(1,1)=aw(1,1)+(dz(j)*rgpt+dz(j+1)*rgp)*indzp
        ELSE
          au(k,k)=au(k,k)+(dr(i)*rlkr(k-1)+dr(i+1)*rlk(k-1))*indrp
          aw(k,k)=aw(k,k)+(rlkt(k-1)*dz(j)+dz(j+1)*rlk(k-1))*indzp
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
      USE grid, ONLY: myijk, myinds
      USE grid, ONLY: r_, t_
      USE tilde_momentum, ONLY: appu, appw
      USE particles_constants, ONLY: rl, inrl
      USE tilde_momentum, ONLY: rug, rwg, rus, rws
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
      REAL*8 :: eps_e, eps_n
!
      imesh = myijk( ip0_jp0_kp0_, ij)
      j  = ( imesh - 1 ) / nr + 1
      i  = MOD( ( imesh - 1 ), nr) + 1
      ijr = myinds(ip1_jp0_kp0_, ij)
      ijt = myinds(ip0_jp1_kp0_, ij)
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
          bw(1)  = rwg(ij)+ dt * indzp *2.D0* ep_n * (p-pt)
        ELSE
          eps_e = (dr(i)*rlkr(k-1) + dr(i+1)*rlk(k-1)) * indrp * inrl(k-1)
          eps_n = (dz(j)*rlkt(k-1) + dz(j+1)*rlk(k-1)) * indzp * inrl(k-1)
          bu(k)  = rus(k-1,ij) + dt * indrp *2.D0* eps_e * (p-pr)
          bw(k)  = rws(k-1,ij) + dt * indzp *2.D0* eps_n * (p-pt)
        ENDIF
!
! ... Implicit terms in the linear system
!
        ks1=k*(k-1)/2
        DO kk=1,k
          ks=ks1+kk
          au(k,kk)=(dr(i)*appu(ks,ijr)+appu(ks,ij)*dr(i+1))*indrp
          au(kk,k)=au(k,kk)
          aw(k,kk)=(dz(j)*appw(ks,ijt)+dz(j+1)*appw(ks,ij))*indzp
          aw(kk,k)=aw(k,kk)
        END DO
        IF(k.EQ.1) THEN
          au(1,1)=au(1,1)+(dr(i)*rgpr+dr(i+1)*rgp)*indrp
          aw(1,1)=aw(1,1)+(dz(j)*rgpt+dz(j+1)*rgp)*indzp
        ELSE
          au(k,k)=au(k,k)+(dr(i)*rlkr(k-1)+dr(i+1)*rlk(k-1))*indrp
          aw(k,k)=aw(k,k)+(rlkt(k-1)*dz(j)+dz(j+1)*rlk(k-1))*indzp
        ENDIF
      END DO
!
      RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
      SUBROUTINE velsk(ug, wg, us, ws, ugm, wgm, usm, wsm, ij)
!
      USE grid, ONLY: fl_l
      USE grid, ONLY: myijk
      USE dimensions
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: ug,wg,us(:),ws(:),ugm,wgm,usm(:),wsm(:)
      INTEGER, INTENT(IN) :: ij
!
      INTEGER :: kk, kp1, k, kj, ki
      INTEGER :: i,j,imesh
      INTEGER :: imj, ijm, ipj, ijp
      REAL*8 :: flt, flr, fll, div, amul, flb
!
      imesh = myijk( ip0_jp0_kp0_, ij)
      j  = ( imesh - 1 ) / nr + 1
      i  = MOD( ( imesh - 1 ), nr) + 1
!
! ... Use Gauss-Jordan method for matrix inversion
!
      imj = myijk( im1_jp0_kp0_, ij)
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
          usm(k-1)=bu1(k)
        END DO

      END IF

      ijm = myijk( ip0_jm1_kp0_, ij)
      flb=fl_l(ijm)
      IF(.NOT.(flb.EQ.2.OR.flb.eq.3.or.flb.eq.5)) THEN

        DO k=2,nphase
          IF(aw1(k,k).LT.rlim) THEN
            DO kk=1,nphase
              aw1(kk,kk)=aw1(kk,kk)+aw1(k,kk)
              aw1(k,kk)=0.D0
              aw1(kk,k)=0.D0
            END DO
            bw1(k)=0.D0
          END IF
        END DO

        DO k=1,nphase
          IF(aw1(k,k).NE.0.D0) THEN
            kp1=k+1
            div=1.D0/aw1(k,k)
            DO kj=kp1,nphase
              aw1(k,kj)=aw1(k,kj)*div
            END DO
            bw1(k)=bw1(k)*div
            aw1(k,k)=0.D0
            DO ki=1,nphase
              amul=aw1(ki,k)
              DO kj=kp1,nphase
                aw1(ki,kj)=aw1(ki,kj)-amul*aw1(k,kj)
              END DO
              bw1(ki)=bw1(ki)-amul*bw1(k)
            END DO
          END IF
        END DO
!
        wgm=bw1(1)
        DO k=2, nphase
          wsm(k-1)=bw1(k)
        END DO
!
      END IF

      ipj = myijk( ip1_jp0_kp0_, ij)
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
          us(k-1)=bu(k)
        END DO
!
      END IF

      ijp = myijk( ip0_jp1_kp0_, ij)
      flt=fl_l(ijp)
      IF(.NOT.(flt.EQ.2.OR.flt.eq.3.or.flt.eq.5)) THEN

        DO k=2,nphase
          IF(aw(k,k).LT.rlim) THEN
            DO kk=1,nphase
              aw(kk,kk)=aw(kk,kk)+aw(k,kk)
              aw(k,kk)=0.D0
              aw(kk,k)=0.D0
            END DO
            bw(k)=0.D0
          END IF
        END DO

        DO k=1,nphase
          IF(aw(k,k).NE.0.D0) THEN 
            kp1=k+1
            div=1.D0/aw(k,k)
            DO kj=kp1,nphase
              aw(k,kj)=aw(k,kj)*div
            END DO
            bw(k)=bw(k)*div
            aw(k,k)=0.D0
            DO ki=1,nphase
              amul=aw(ki,k)
              DO kj=kp1,nphase
                aw(ki,kj)=aw(ki,kj)-amul*aw(k,kj)
              END DO 
              bw(ki)=bw(ki)-amul*bw(k)
            END DO 
          END IF
        END DO
!
        wg=bw(1)
        DO k=2,nphase
          ws(k-1)=bw(k)
        END DO
!
      END IF
!
      RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
      SUBROUTINE velsk2(ug, wg, us, ws, ij)
!
      USE grid, ONLY: fl_l
      USE grid, ONLY: myijk
      USE dimensions
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: ug, wg, us(:), ws(:)
      INTEGER, INTENT(IN) :: ij
      INTEGER :: ipj, ijp
!
      INTEGER :: kk, kp1, k, kj, ki
      REAL*8 :: flt, flr, fll, div, amul, flb

      ipj = myijk( ip1_jp0_kp0_, ij)
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
          us(k-1)=bu(k)
        END DO
!
      END IF

      ijp = myijk( ip0_jp1_kp0_, ij)
      flt=fl_l(ijp)
      IF(.NOT.(flt.EQ.2.OR.flt.eq.3.or.flt.eq.5)) THEN

        DO k=2,nphase
          IF(aw(k,k).LT.rlim) THEN
            DO kk=1,nphase
              aw(kk,kk)=aw(kk,kk)+aw(k,kk)
              aw(k,kk)=0.D0
              aw(kk,k)=0.D0
            END DO
            bw(k)=0.D0
          END IF
        END DO

        DO k=1,nphase
          IF(aw(k,k).NE.0.D0) THEN 
            kp1=k+1
            div=1.D0/aw(k,k)
            DO kj=kp1,nphase
              aw(k,kj)=aw(k,kj)*div
            END DO
            bw(k)=bw(k)*div
            aw(k,k)=0.D0
            DO ki=1,nphase
              amul=aw(ki,k)
              DO kj=kp1,nphase
                aw(ki,kj)=aw(ki,kj)-amul*aw(k,kj)
              END DO 
              bw(ki)=bw(ki)-amul*bw(k)
            END DO 
          END IF
        END DO
!
        wg=bw(1)
        DO k=2,nphase
          ws(k-1)=bw(k)
        END DO
!
      END IF
!
      RETURN
      END SUBROUTINE

!------------------------------------------------------------------------
      END MODULE
!------------------------------------------------------------------------
