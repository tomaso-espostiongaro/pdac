!--------------------------------------------------------------------
      MODULE phases_matrix
!--------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

      REAL*8, PRIVATE, DIMENSION(:),   ALLOCATABLE ::  bu, bv, bw
      REAL*8, PRIVATE, DIMENSION(:,:), ALLOCATABLE ::  au, av, aw
      REAL*8, PRIVATE, DIMENSION(:),   ALLOCATABLE ::  bu1, bv1, bw1
      REAL*8, PRIVATE, DIMENSION(:,:), ALLOCATABLE ::  au1, av1, aw1

      REAL*8 :: rlim
!--------------------------------------------------------------------
      CONTAINS
!--------------------------------------------------------------------
      SUBROUTINE bounds_matrix
      USE dimensions
      IMPLICIT NONE
!
      ALLOCATE(bu1(nphase), bv1(nphase), bw1(nphase))
      ALLOCATE(bu(nphase), bv(nphase), bw(nphase))
      ALLOCATE(au1(nphase,nphase), av1(nphase,nphase), aw1(nphase,nphase))
      ALLOCATE(au(nphase,nphase), av(nphase,nphase), aw(nphase,nphase))
!
      au = 0.D0;   av = 0.D0;   aw = 0.D0
      au1 = 0.D0;  av1 = 0.D0;  aw1 = 0.D0
      bu = 0.D0;   bv = 0.D0;   bw = 0.D0
      bu1 = 0.D0;  bv1 = 0.D0;  bw1 = 0.D0
!
      RETURN
      END SUBROUTINE
!--------------------------------------------------------------------
      SUBROUTINE mats(ijk)
!
! ... Computes matrix elements to solve momentum-balance 
! ... linear system of coupled equations, in current cell and in
! ... west, south and bottom cells 
!
      USE dimensions
      USE gas_solid_density, ONLY: rgp, rlk
      USE grid, ONLY: dx, dy, dz
      USE grid, ONLY: myijk
      USE indijk_module
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p, ep
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes, ONLY: ijke, ijkn, ijkt, ijkw, ijks, ijkb
      USE tilde_momentum, ONLY: rug, rvg, rwg, rus, rvs, rws
      USE tilde_momentum, ONLY: appu, appv,  appw
      USE time_parameters, ONLY: dt

      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ijk
!
      INTEGER :: ls, ll, ls1, l
      INTEGER :: i,j,k,imesh,is
      REAL*8 :: dxp, dyp, dzp, dxm, dym, dzm
      REAL*8 :: indxp, indyp, indzp, indxm, indym, indzm
      REAL*8 :: ep_w, ep_s, ep_b, ep_e, ep_n, ep_t
      REAL*8 :: eps_w, eps_s, eps_b, eps_e, eps_n, eps_t
      REAL*8 :: dxi, dxim1, dyj, dyjm1, dzk, dzkm1
      REAL*8 :: dxip1, dyjp1, dzkp1
      REAL*8 :: pijk
!
      imesh = myijk( ip0_jp0_kp0_, ijk)
      i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
      j = MOD( imesh - 1, nx*ny ) / nx + 1
      k = ( imesh - 1 ) / ( nx*ny ) + 1
!
      dxm=dx(i)+dx(i-1)
      dym=dy(j)+dy(j-1)
      dzm=dz(k)+dz(k-1)
      dxi=dx(i)
      dxim1=dx(i-1)
      dxip1=dx(i+1)
      dyj=dy(j)
      dyjm1=dy(j-1)
      dyjp1=dy(j+1)
      dzk=dz(k)
      dzkm1=dz(k-1)
      dzkp1=dz(k+1)
      pijk = p(ijk)
!
      indxm=1.D0/dxm
      indym=1.D0/dym
      indzm=1.D0/dzm
!
! ... Backward ...
!
      ep_w = (dxi*ep(ijkw) + dxim1*ep(ijk)) * indxm
      ep_s = (dyj*ep(ijks) + dyjm1*ep(ijk)) * indym
      ep_b = (dzk*ep(ijkb) + dzkm1*ep(ijk)) * indzm
      bu1(1) = rug(imjk)+ dt * indxm *2.D0* ep_w * (p(ijkw)-pijk)
      bv1(1) = rvg(ijmk)+ dt * indym *2.D0* ep_s * (p(ijks)-pijk)
      bw1(1) = rwg(ijkm)+ dt * indzm *2.D0* ep_b * (p(ijkb)-pijk)
!
      au1(1,1)=(dxi*appu(1,ijkw)+dxim1*appu(1,ijk))*indxm
      av1(l,1)=(dyj*appv(1,ijks)+dyjm1*appv(1,ijk))*indym
      aw1(1,1)=(dzk*appw(1,ijkb)+dzkm1*appw(1,ijk))*indzm
      au1(1,1)=au1(1,1)+(dxi*rgp(ijkw)+dxim1*rgp(ijk))*indxm
      av1(1,1)=av1(1,1)+(dyj*rgp(ijks)+dyjm1*rgp(ijk))*indym
      aw1(1,1)=aw1(1,1)+(dzk*rgp(ijkb)+dzkm1*rgp(ijk))*indzm

      DO l=2,nphase
!
! ... Explicit terms in the linear system
!
          eps_w = (dxi*rlk(ijkw,l-1) + dxim1*rlk(ijk,l-1)) * indxm * inrl(l-1)
          eps_s = (dyj*rlk(ijks,l-1) + dyjm1*rlk(ijk,l-1)) * indym * inrl(l-1)
          eps_b = (dzk*rlk(ijkb,l-1) + dzkm1*rlk(ijk,l-1)) * indzm * inrl(l-1)
          bu1(l) = rus(l-1,imjk) + dt * indxm *2.D0* eps_w * (p(ijkw)-pijk)
          bv1(l) = rws(l-1,ijmk) + dt * indym *2.D0* eps_s * (p(ijks)-pijk)
          bw1(l) = rws(l-1,ijkm) + dt * indzm *2.D0* eps_b * (p(ijkb)-pijk)
!
! ... Implicit terms in the linear system
!
! ... number of elements of appu, appv, appw, already used in previous cycles
! ... (= size of the upper triangular matrix of side l)
        ls1=l*(l-1)/2
!
! ... Off-diagonal elements
        DO ll=1,l
          ls=ls1+ll
          au1(l,ll)=(dxi*appu(ls,ijkw)+dxim1*appu(ls,ijk))*indxm
          au1(ll,l)=au1(l,ll)
          av1(l,ll)=(dyj*appv(ls,ijks)+dyjm1*appv(ls,ijk))*indym
          av1(ll,l)=av1(l,ll)
          aw1(l,ll)=(dzk*appw(ls,ijkb)+dzkm1*appw(ls,ijk))*indzm
          aw1(ll,l)=aw1(l,ll)
        END DO
!
! ... Diagonal elements

          au1(l,l)=au1(l,l)+(dxi*rlk(ijkw,l-1)+dxim1*rlk(ijk,l-1))*indxm
          av1(l,l)=av1(l,l)+(dyj*rlk(ijks,l-1)+dyjm1*rlk(ijk,l-1))*indym
          aw1(l,l)=aw1(l,l)+(dzk*rlk(ijkb,l-1)+dzkm1*rlk(ijk,l-1))*indzm

      END DO
!
! ... Forward ...
!
      dxp=dxi+dx(i+1)
      dyp=dyj+dy(j+1)
      dzp=dzk+dz(k+1)
!
      indxp=1.D0/dxp
      indyp=1.D0/dyp
      indzp=1.D0/dzp
!
      ep_e = (dxi*ep(ijke) + dxip1*ep(ijk)) * indxp
      ep_n = (dyj*ep(ijkn) + dyjp1*ep(ijk)) * indyp
      ep_t = (dzk*ep(ijkt) + dzkp1*ep(ijk)) * indzp
      bu(1)  = rug(ijk)+ dt * indxp *2.D0* ep_e * (pijk-p(ijke))
      bv(1)  = rvg(ijk)+ dt * indyp *2.D0* ep_n * (pijk-p(ijkn))
      bw(1)  = rwg(ijk)+ dt * indzp *2.D0* ep_t * (pijk-p(ijkt))
!
      au(1,1)=(dxi*appu(1,ijke)+dxip1*appu(1,ijk))*indxp
      av(1,1)=(dyj*appv(1,ijkn)+dyjp1*appv(1,ijk))*indyp
      aw(1,1)=(dzk*appw(1,ijkt)+dzkp1*appw(1,ijk))*indzp
      au(1,1)=au(1,1)+(dxi*rgp(ijke)+dxip1*rgp(ijk))*indxp
      av(1,1)=av(1,1)+(dyj*rgp(ijkn)+dyjp1*rgp(ijk))*indyp
      aw(1,1)=aw(1,1)+(dzk*rgp(ijkt)+dzkp1*rgp(ijk))*indzp

      DO l=2,nphase
!
! ... Explicit terms in the linear system
!
        eps_e = (dxi*rlk(ijke,l-1) + dxip1*rlk(ijk,l-1)) * indxp * inrl(l-1)
        eps_n = (dyj*rlk(ijkn,l-1) + dyjp1*rlk(ijk,l-1)) * indyp * inrl(l-1)
        eps_t = (dzk*rlk(ijkt,l-1) + dzkp1*rlk(ijk,l-1)) * indzp * inrl(l-1)
        bu(l)  = rus(l-1,ijk) + dt * indxp *2.D0* eps_e * (pijk-p(ijke))
        bv(l)  = rvs(l-1,ijk) + dt * indyp *2.D0* eps_n * (pijk-p(ijkn))
        bw(l)  = rws(l-1,ijk) + dt * indzp *2.D0* eps_t * (pijk-p(ijkt))
!
! ... Implicit terms in the linear system
!
! ... number of elements of appu, appv, appw already used in previous cycles
! ... (= size of the upper triangular matrix of side l)
        ls1=l*(l-1)/2
!
! ... Off-Diagonal elements
        DO ll=1,l
          ls=ls1+ll
          au(l,ll)=(dxi*appu(ls,ijke)+dxip1*appu(ls,ijk))*indxp
          au(ll,l)=au(l,ll)
          av(l,ll)=(dyj*appv(ls,ijkn)+dyjp1*appv(ls,ijk))*indyp
          av(ll,l)=av(l,ll)
          aw(l,ll)=(dzk*appw(ls,ijkt)+dzkp1*appw(ls,ijk))*indzp
          aw(ll,l)=aw(l,ll)
        END DO
!
! ... Diagonal elements
        au(l,l)=au(l,l)+(dxi*rlk(ijke,l-1)+dxip1*rlk(ijk,l-1))*indxp
        av(l,l)=av(l,l)+(dyj*rlk(ijkn,l-1)+dyjp1*rlk(ijk,l-1))*indyp
        aw(l,l)=aw(l,l)+(dzk*rlk(ijkt,l-1)+dzkp1*rlk(ijk,l-1))*indzp

      END DO
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE mats2(ijk)
!
! ... Computes matrix elements to solve momentum-balance 
! ... linear system of coupled equations in current cell
!
      USE dimensions
      USE gas_solid_density, ONLY: rgp, rlk
      USE grid, ONLY: dx, dy, dz
      USE grid, ONLY: myijk
      USE indijk_module
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p, ep
      USE set_indexes, ONLY: ijke, ijkn, ijkt
      USE tilde_momentum, ONLY: rug, rvg, rwg, rus, rvs, rws
      USE tilde_momentum, ONLY: appu, appv, appw
      USE time_parameters, ONLY: dt

      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ijk
!
      INTEGER :: ls, ll, ls1, l
      INTEGER :: i,j,k,imesh,is
      REAL*8 :: dxp, dyp, dzp, indxp, indyp, indzp
      REAL*8 :: ep_e, ep_n, ep_t
      REAL*8 :: eps_e, eps_n, eps_t
      REAL*8 :: dxi, dxim1, dyj, dyjm1, dzk, dzkm1
      REAL*8 :: dxip1, dyjp1, dzkp1
      REAL*8 :: pijk
!
      imesh = myijk( ip0_jp0_kp0_, ijk)
      i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
      j = MOD( imesh - 1, nx*ny ) / nx + 1
      k = ( imesh - 1 ) / ( nx*ny ) + 1
!
      dxi=dx(i)
      dxip1=dx(i+1)
      dyj=dy(j)
      dyjp1=dy(j+1)
      dzk=dz(k)
      dzkp1=dz(k+1)
      pijk = p(ijk)
!
! ... Forward ...
!
      dxp=dxi+dxip1
      dyp=dyj+dyjp1
      dzp=dzk+dzkp1
!
      indxp=1.D0/dxp
      indyp=1.D0/dyp
      indzp=1.D0/dzp
!
      ep_e = (dxi*ep(ijke) + dxip1*ep(ijk)) * indxp
      ep_n = (dyj*ep(ijkn) + dyjp1*ep(ijk)) * indyp
      ep_t = (dzk*ep(ijkt) + dzkp1*ep(ijk)) * indzp
      bu(1)  = rug(ijk)+ dt * indxp *2.D0* ep_e * (pijk-p(ijke))
      bv(1)  = rvg(ijk)+ dt * indyp *2.D0* ep_n * (pijk-p(ijkn))
      bw(1)  = rwg(ijk)+ dt * indzp *2.D0* ep_t * (pijk-p(ijkt))
!
      au(1,1)=(dxi*appu(1,ijke)+dxip1*appu(1,ijk))*indxp
      av(1,1)=(dyj*appv(1,ijkn)+dyjp1*appv(1,ijk))*indyp
      aw(1,1)=(dzk*appw(1,ijkt)+dzkp1*appw(1,ijk))*indzp
      au(1,1)=au(1,1)+(dxi*rgp(ijke)+dxip1*rgp(ijk))*indxp
      av(1,1)=av(1,1)+(dyj*rgp(ijkn)+dyjp1*rgp(ijk))*indyp
      aw(1,1)=aw(1,1)+(dzk*rgp(ijkt)+dzkp1*rgp(ijk))*indzp

      DO l=2,nphase
!
! ... Explicit terms in the linear system
!
        eps_e = (dxi*rlk(ijke,l-1) + dxip1*rlk(ijk,l-1)) * indxp * inrl(l-1)
        eps_n = (dyj*rlk(ijkn,l-1) + dyjp1*rlk(ijk,l-1)) * indyp * inrl(l-1)
        eps_t = (dzk*rlk(ijkt,l-1) + dzkp1*rlk(ijk,l-1)) * indzp * inrl(l-1)
        bu(l)  = rus(l-1,ijk) + dt * indxp *2.D0* eps_e * (pijk-p(ijke))
        bv(l)  = rvs(l-1,ijk) + dt * indyp *2.D0* eps_n * (pijk-p(ijkn))
        bw(l)  = rws(l-1,ijk) + dt * indzp *2.D0* eps_t * (pijk-p(ijkt))
!
! ... Implicit terms in the linear system
!
! ... number of elements of appu, appv, appw already used in previous cycles
! ... (= size of the upper triangular matrix of side l)
        ls1=l*(l-1)/2
!
! ... Off-Diagonal elements
        DO ll=1,l
          ls=ls1+ll
          au(l,ll)=(dxi*appu(ls,ijke)+dxip1*appu(ls,ijk))*indxp
          au(ll,l)=au(l,ll)
          av(l,ll)=(dyj*appv(ls,ijkn)+dyjp1*appv(ls,ijk))*indyp
          av(ll,l)=av(l,ll)
          aw(l,ll)=(dzk*appw(ls,ijkt)+dzkp1*appw(ls,ijk))*indzp
          aw(ll,l)=aw(l,ll)
        END DO
!
! ... Diagonal elements
        au(l,l)=au(l,l)+(dxi*rlk(ijke,l-1)+dxip1*rlk(ijk,l-1))*indxp
        av(l,l)=av(l,l)+(dyj*rlk(ijkn,l-1)+dyjp1*rlk(ijk,l-1))*indyp
        aw(l,l)=aw(l,l)+(dzk*rlk(ijkt,l-1)+dzkp1*rlk(ijk,l-1))*indzp

      END DO
!
      RETURN
      END SUBROUTINE mats2
!----------------------------------------------------------------------
      SUBROUTINE velsk(ug, vg, wg, us, vs, ws, ugm, vgm, wgm, usm, vsm, wsm)
!
      USE dimensions, ONLY: nphase
      USE grid, ONLY: fl_l
      USE set_indexes, ONLY: imjk, ijmk, ijkm, ipjk, ijpk, ijkp
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: ug,vg,wg,us(:),vs(:),ws(:)
      REAL*8, INTENT(OUT) :: ugm,vgm,wgm,usm(:),vsm(:),wsm(:)
!
      INTEGER :: ll, lp1, l, lj, li
      INTEGER :: flw, fls, flb, fle, fln, flt
      REAL*8 :: div, amul
!
! ... Use Gauss-Jordan method for matrix inversion
!
      flw=fl_l(imjk)
      IF (flw == 1 .OR. flw == 4) THEN

        DO l=2,nphase
          IF(au1(l,l) < rlim) THEN
            DO ll=1,nphase
              au1(ll,ll)=au1(ll,ll)+au1(l,ll)
              au1(l,ll)=0.D0
              au1(ll,l)=0.D0
            END DO
            bu1(l)=0.D0
          END IF
        END DO

        DO l=1,nphase
          IF(au1(l,l) /= 0.D0) THEN 
            lp1=l+1
            div=1.D0/au1(l,l)
            DO lj=lp1,nphase
              au1(l,lj)=au1(l,lj)*div
            END DO
            bu1(l)=bu1(l)*div
            au1(l,l)=0.D0
            DO li=1,nphase
              amul=au1(li,l)
              DO lj=lp1,nphase
                au1(li,lj)=au1(li,lj)-amul*au1(l,lj)
              END DO
              bu1(li)=bu1(li)-amul*bu1(l)
            END DO
          END IF
        END DO
!
        ugm=bu1(1)
        DO l=2,nphase
          usm(l-1)=bu1(l)
        END DO

      END IF

      fls=fl_l(ijmk)
      IF (fls == 1 .OR. fls == 4) THEN

        DO l=2,nphase
          IF(av1(l,l) < rlim) THEN
            DO ll=1,nphase
              av1(ll,ll)=av1(ll,ll)+av1(l,ll)
              av1(l,ll)=0.D0
              av1(ll,l)=0.D0
            END DO
            bv1(l)=0.D0
          END IF
        END DO

        DO l=1,nphase
          IF(av1(l,l) /= 0.D0) THEN
            lp1=l+1
            div=1.D0/av1(l,l)
            DO lj=lp1,nphase
              av1(l,lj)=av1(l,lj)*div
            END DO
            bv1(l)=bv1(l)*div
            av1(l,l)=0.D0
            DO li=1,nphase
              amul=av1(li,l)
              DO lj=lp1,nphase
                av1(li,lj)=av1(li,lj)-amul*av1(l,lj)
              END DO
              bv1(li)=bv1(li)-amul*bv1(l)
            END DO
          END IF
        END DO
!
        vgm=bv1(1)
        DO l=2, nphase
          vsm(l-1)=bv1(l)
        END DO
!
      END IF

      flb=fl_l(ijkm)
      IF (flb == 1 .OR. flb == 4) THEN

        DO l=2,nphase
          IF(aw1(l,l) < rlim) THEN
            DO ll=1,nphase
              aw1(ll,ll)=aw1(ll,ll)+aw1(l,ll)
              aw1(l,ll)=0.D0
              aw1(ll,l)=0.D0
            END DO
            bw1(l)=0.D0
          END IF
        END DO

        DO l=1,nphase
          IF(aw1(l,l) /= 0.D0) THEN
            lp1=l+1
            div=1.D0/aw1(l,l)
            DO lj=lp1,nphase
              aw1(l,lj)=aw1(l,lj)*div
            END DO
            bw1(l)=bw1(l)*div
            aw1(l,l)=0.D0
            DO li=1,nphase
              amul=aw1(li,l)
              DO lj=lp1,nphase
                aw1(li,lj)=aw1(li,lj)-amul*aw1(l,lj)
              END DO
              bw1(li)=bw1(li)-amul*bw1(l)
            END DO
          END IF
        END DO
!
        wgm=bw1(1)
        DO l=2, nphase
          wsm(l-1)=bw1(l)
        END DO
!
      END IF

      fle=fl_l(ipjk)
      IF (fle == 1 .OR. fle == 4) THEN

        DO l=2,nphase
          IF(au(l,l) < rlim) THEN
            DO ll=1,nphase
              au(ll,ll)=au(ll,ll)+au(l,ll)
              au(l,ll)=0.D0
              au(ll,l)=0.D0
            END DO
            bu(l)=0.D0
          END IF
        END DO

        DO l=1,nphase
          IF(au(l,l) /= 0.D0) THEN
            lp1=l+1
            div=1.D0/au(l,l)
            DO lj=lp1,nphase
              au(l,lj)=au(l,lj)*div
            END DO
            bu(l)=bu(l)*div
            au(l,l)=0.D0
            DO li=1,nphase
              amul=au(li,l)
              DO lj=lp1,nphase
                au(li,lj)=au(li,lj)-amul*au(l,lj)
              END DO
              bu(li)=bu(li)-amul*bu(l)
            END DO
          END IF
        END DO
!
        ug=bu(1)
        DO l=2,nphase
          us(l-1)=bu(l)
        END DO
!
      END IF

      fln=fl_l(ijpk)
      IF (fln == 1 .OR. fln == 4) THEN
        DO l=2,nphase
          IF(av(l,l) < rlim) THEN
            DO ll=1,nphase
              av(ll,ll)=av(ll,ll)+av(l,ll)
              av(l,ll)=0.D0
              av(ll,l)=0.D0
            END DO
            bv(l)=0.D0
          END IF
        END DO

        DO l=1,nphase
          IF(av(l,l) /= 0.D0) THEN 
            lp1=l+1
            div=1.D0/av(l,l)
            DO lj=lp1,nphase
              av(l,lj)=av(l,lj)*div
            END DO
            bv(l)=bv(l)*div
            av(l,l)=0.D0
            DO li=1,nphase
              amul=av(li,l)
              DO lj=lp1,nphase
                av(li,lj)=av(li,lj)-amul*av(l,lj)
              END DO 
              bv(li)=bv(li)-amul*bv(l)
            END DO 
          END IF
        END DO
!
        vg=bv(1)
        DO l=2,nphase
          vs(l-1)=bv(l)
        END DO
!
      END IF
!

      flt=fl_l(ijkp)
      IF (flt == 1 .OR. flt == 4) THEN

        DO l=2,nphase
          IF(aw(l,l) < rlim) THEN
            DO ll=1,nphase
              aw(ll,ll)=aw(ll,ll)+aw(l,ll)
              aw(l,ll)=0.D0
              aw(ll,l)=0.D0
            END DO
            bw(l)=0.D0
          END IF
        END DO

        DO l=1,nphase
          IF(aw(l,l) /= 0.D0) THEN 
            lp1=l+1
            div=1.D0/aw(l,l)
            DO lj=lp1,nphase
              aw(l,lj)=aw(l,lj)*div
            END DO
            bw(l)=bw(l)*div
            aw(l,l)=0.D0
            DO li=1,nphase
              amul=aw(li,l)
              DO lj=lp1,nphase
                aw(li,lj)=aw(li,lj)-amul*aw(l,lj)
              END DO 
              bw(li)=bw(li)-amul*bw(l)
            END DO 
          END IF
        END DO
!
        wg=bw(1)
        DO l=2,nphase
          ws(l-1)=bw(l)
        END DO
!
      END IF
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE velsk2(ijk)
!
      USE dimensions, ONLY: nphase
      USE grid, ONLY: fl_l
      USE set_indexes, ONLY: ipjk, ijpk, ijkp
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ijk
!
      INTEGER :: ll, lp1, l, lj, li
      INTEGER :: fle, fln, flt
      REAL*8 :: div, amul
!
! ... Use Gauss-Jordan method for matrix inversion
!
      fle=fl_l(ipjk)
      IF (fle == 1 .OR. fle == 4) THEN

        DO l=2,nphase
          IF(au(l,l) < rlim) THEN
            DO ll=1,nphase
              au(ll,ll)=au(ll,ll)+au(l,ll)
              au(l,ll)=0.D0
              au(ll,l)=0.D0
            END DO
            bu(l)=0.D0
          END IF
        END DO

        DO l=1,nphase
          IF(au(l,l) /= 0.D0) THEN
            lp1=l+1
            div=1.D0/au(l,l)
            DO lj=lp1,nphase
              au(l,lj)=au(l,lj)*div
            END DO
            bu(l)=bu(l)*div
            au(l,l)=0.D0
            DO li=1,nphase
              amul=au(li,l)
              DO lj=lp1,nphase
                au(li,lj)=au(li,lj)-amul*au(l,lj)
              END DO
              bu(li)=bu(li)-amul*bu(l)
            END DO
          END IF
        END DO
!
        ug(ijk)=bu(1)
        DO l=2,nphase
          us(ijk,l-1)=bu(l)
        END DO
!
      END IF

      fln=fl_l(ijpk)
      IF (fln == 1 .OR. fln == 4) THEN
        DO l=2,nphase
          IF(av(l,l) < rlim) THEN
            DO ll=1,nphase
              av(ll,ll)=av(ll,ll)+av(l,ll)
              av(l,ll)=0.D0
              av(ll,l)=0.D0
            END DO
            bv(l)=0.D0
          END IF
        END DO

        DO l=1,nphase
          IF(av(l,l) /= 0.D0) THEN 
            lp1=l+1
            div=1.D0/av(l,l)
            DO lj=lp1,nphase
              av(l,lj)=av(l,lj)*div
            END DO
            bv(l)=bv(l)*div
            av(l,l)=0.D0
            DO li=1,nphase
              amul=av(li,l)
              DO lj=lp1,nphase
                av(li,lj)=av(li,lj)-amul*av(l,lj)
              END DO 
              bv(li)=bv(li)-amul*bv(l)
            END DO 
          END IF
        END DO
!
        vg(ijk)=bv(1)
        DO l=2,nphase
          vs(ijk,l-1)=bv(l)
        END DO
!
      END IF
!

      flt=fl_l(ijkp)
      IF (flt == 1 .OR. flt == 4) THEN

        DO l=2,nphase
          IF(aw(l,l) < rlim) THEN
            DO ll=1,nphase
              aw(ll,ll)=aw(ll,ll)+aw(l,ll)
              aw(l,ll)=0.D0
              aw(ll,l)=0.D0
            END DO
            bw(l)=0.D0
          END IF
        END DO

        DO l=1,nphase
          IF(aw(l,l) /= 0.D0) THEN 
            lp1=l+1
            div=1.D0/aw(l,l)
            DO lj=lp1,nphase
              aw(l,lj)=aw(l,lj)*div
            END DO
            bw(l)=bw(l)*div
            aw(l,l)=0.D0
            DO li=1,nphase
              amul=aw(li,l)
              DO lj=lp1,nphase
                aw(li,lj)=aw(li,lj)-amul*aw(l,lj)
              END DO 
              bw(li)=bw(li)-amul*bw(l)
            END DO 
          END IF
        END DO
!
        wg(ijk)=bw(1)
        DO l=2,nphase
          ws(ijk,l-1)=bw(l)
        END DO
!
      END IF
!
      RETURN
      END SUBROUTINE
!------------------------------------------------------------------------
      END MODULE
!------------------------------------------------------------------------
