!--------------------------------------------------------------------
      MODULE phases_matrix
!--------------------------------------------------------------------
      USE particles_constants, ONLY: plim
      IMPLICIT NONE
      SAVE
!
! ... Interphase momentum matrix elements
!
      REAL*8, PRIVATE, DIMENSION(:),   ALLOCATABLE ::  bu, bv, bw
      REAL*8, PRIVATE, DIMENSION(:,:), ALLOCATABLE ::  au, av, aw
      REAL*8, PRIVATE, DIMENSION(:),   ALLOCATABLE ::  bu1, bv1, bw1
      REAL*8, PRIVATE, DIMENSION(:,:), ALLOCATABLE ::  au1, av1, aw1

      REAL*8 :: rlim

!--------------------------------------------------------------------
      CONTAINS
!--------------------------------------------------------------------
      SUBROUTINE allocate_matrix
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
      END SUBROUTINE allocate_matrix

!--------------------------------------------------------------------
      SUBROUTINE assemble_all_matrix(ijk)
!
! ... Computes matrix elements to solve momentum-balance 
! ... linear system of coupled equations on all cell faces
!
      USE atmospheric_conditions, ONLY: gravz, gravx, gravy
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions
      USE domain_mapping, ONLY: myijk, meshinds
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE grid, ONLY: dx, dy, dz
      USE indijk_module
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p, ep, pmodel, pmodel
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes, ONLY: ijkw, ijks, ijkb
      USE tilde_momentum, ONLY: rug, rvg, rwg, rus, rvs, rws
      USE tilde_momentum, ONLY: appu, appv, appw
      USE time_parameters, ONLY: dt, time, alpha, alphagrav

      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ijk
!
      INTEGER :: ls, ll, ls1, l
      INTEGER :: i,j,k,imesh,is
      REAL*8 :: dxm, dym, dzm
      REAL*8 :: indxm, indym, indzm
      REAL*8 :: ep_w, ep_s, ep_b
      REAL*8 :: eps_w, eps_s, eps_b
      REAL*8 :: dxi, dxim1, dyj, dyjm1, dzk, dzkm1
      REAL*8 :: pijk 
      REAL*8 :: epijk, rgpijk, rlklm1 
!
      CALL meshinds(ijk,imesh,i,j,k)
!
      dxi   = dx(i)
      dxm   = dx(i)+dx(i-1)
      dxim1 = dx(i-1)
      dyj   = dy(j)
      dym   = dy(j)+dy(j-1)
      dyjm1 = dy(j-1)
      dzk   = dz(k)
      dzm   = dz(k)+dz(k-1)
      dzkm1 = dz(k-1)
      pijk  = p(ijk)
      epijk = ep(ijk)
      rgpijk= rgp(ijk)
!
      indxm=1.D0/dxm
      indzm=1.D0/dzm
      indym=1.D0/dym
!
      IF ( pmodel == 1 ) THEN
         ep_w = (dxi*ep(ijkw) + dxim1*epijk) * indxm
         ep_s = (dyj*ep(ijks) + dyjm1*epijk) * indym
         ep_b = (dzk*ep(ijkb) + dzkm1*epijk) * indzm
      ELSE IF ( pmodel == 2 ) THEN
         ep_w = 1.D0
         ep_s = 1.D0
         ep_b = 1.D0
      END IF
!
      bu1(1) = rug(imjk) + alpha * dt * indxm *2.D0* ep_w * (p(ijkw)-pijk)
      bw1(1) = rwg(ijkm) + alpha * dt * indzm *2.D0* ep_b * (p(ijkb)-pijk)
!
      bu1(1) = bu1(1) + dt*alphagrav*indxm*gravx*(dxi*rgp(ijkw) + dxim1*rgp(ijk))
      bw1(1) = bw1(1) + dt*alphagrav*indzm*gravz*(dzk*rgp(ijkb) + dzkm1*rgp(ijk))
!
      au1(1,1)=alpha*(dxi*appu(ijkw,1)+dxim1*appu(ijk,1))*indxm
      au1(1,1)=au1(1,1)+(dxi*rgp(ijkw)+dxim1*rgpijk)*indxm
      aw1(1,1)=alpha*(dzk*appw(ijkb,1)+dzkm1*appw(ijk,1))*indzm
      aw1(1,1)=aw1(1,1)+(dzk*rgp(ijkb)+dzkm1*rgpijk)*indzm
!
      IF (job_type == JOB_TYPE_3D) THEN

        bv1(1) = rvg(ijmk) + alpha * dt * indym *2.D0* ep_s * (p(ijks)-pijk)
        bv1(1) = bv1(1) + dt*alphagrav*indym*gravy*(dyj*rgp(ijks) + dyjm1*rgp(ijk))
        av1(1,1)=alpha*(dyj*appv(ijks,1)+dyjm1*appv(ijk,1))*indym
        av1(1,1)=av1(1,1) + (dyj*rgp(ijks)+dyjm1*rgpijk)*indym

      END IF
!
      DO l = 2, nphase
!
! ... Explicit terms in the linear system
!
        rlklm1  = rlk(ijk,l-1)

        IF ( pmodel == 1 ) THEN
          eps_w = (dxi*rlk(ijkw,l-1) + dxim1*rlklm1) * indxm * inrl(l-1)
          eps_s = (dyj*rlk(ijks,l-1) + dyjm1*rlklm1) * indym * inrl(l-1)
          eps_b = (dzk*rlk(ijkb,l-1) + dzkm1*rlklm1) * indzm * inrl(l-1)
        ELSE IF ( pmodel == 2 ) THEN
          eps_w = 0.D0
          eps_s = 0.D0
          eps_b = 0.D0
        END IF
!
        bu1(l) = rus(imjk,l-1) + alpha * dt * indxm *2.D0* eps_w * (p(ijkw)-pijk)
        bw1(l) = rws(ijkm,l-1) + alpha * dt * indzm *2.D0* eps_b * (p(ijkb)-pijk)

        bu1(l) = bu1(l) + dt*alphagrav*indxm*gravx*(dxi*rlk(ijkw,l-1) + dxim1*rlk(ijk,l-1))
        bw1(l) = bw1(l) + dt*alphagrav*indzm*gravz*(dzk*rlk(ijkb,l-1) + dzkm1*rlk(ijk,l-1))
!
        IF (job_type == JOB_TYPE_3D) THEN
          bv1(l) = rvs(ijmk,l-1) + alpha * dt * indym *2.D0* eps_s * (p(ijks)-pijk)
          bv1(l) = bv1(l) + dt*alphagrav*indym*gravy*(dyj*rlk(ijks,l-1)+dyjm1*rlk(ijk,l-1))
        END IF
!
! ... Implicit terms in the linear system
!
! ... number of elements of appu, appv, appw, already used in previous cycles
! ... (= size of the upper triangular matrix of side l)
        ls1 = l * (l-1) / 2
!
! ... Off-diagonal elements

        DO ll = 1, l
          ls = ls1 + ll
          au1(ll,l)=alpha*(dxi*appu(ijkw,ls)+dxim1*appu(ijk,ls))*indxm
          au1(l,ll)=au1(ll,l)
          aw1(ll,l)=alpha*(dzk*appw(ijkb,ls)+dzkm1*appw(ijk,ls))*indzm
          aw1(l,ll)=aw1(ll,l)

          IF (job_type == JOB_TYPE_3D) THEN
            av1(ll,l)=alpha*(dyj*appv(ijks,ls)+dyjm1*appv(ijk,ls))*indym
            av1(l,ll)=av1(ll,l)
          END IF
        END DO
!
! ... Diagonal elements

        au1(l,l)=au1(l,l)+(dxi*rlk(ijkw,l-1)+dxim1*rlklm1)*indxm
        aw1(l,l)=aw1(l,l)+(dzk*rlk(ijkb,l-1)+dzkm1*rlklm1)*indzm

        IF (job_type == JOB_TYPE_3D) av1(l,l)=av1(l,l)+(dyj*rlk(ijks,l-1)+dyjm1*rlklm1)*indym

      END DO
!
      CALL assemble_matrix(ijk)
!
      RETURN
      END SUBROUTINE assemble_all_matrix
!--------------------------------------------------------------------
      SUBROUTINE assemble_matrix(ijk)
!
! ... Computes matrix elements to solve momentum-balance 
! ... linear system of coupled equations 
! ... only on East, North, and Top faces of the cell
!
      USE atmospheric_conditions, ONLY: gravz, gravx, gravy
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions
      USE domain_mapping, ONLY: myijk, meshinds
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: tg, ts
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE grid, ONLY: dx, dy, dz
      USE indijk_module
      USE io_files, ONLY: testunit
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p, ep, pmodel
      USE set_indexes, ONLY: ijke, ijkn, ijkt
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE tilde_momentum, ONLY: rug, rvg, rwg, rus, rvs, rws
      USE tilde_momentum, ONLY: appu, appv, appw
      USE time_parameters, ONLY: dt, time, alpha, alphagrav

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
      REAL*8 :: pijk, epijk
!
      CALL meshinds(ijk,imesh,i,j,k)

      dxi   = dx(i)
      dxip1 = dx(i+1)
      dyj   = dy(j)
      dyjp1 = dy(j+1)
      dzk   = dz(k)
      dzkp1 = dz(k+1)
      pijk  = p(ijk)
      epijk = ep(ijk)
!
      dxp=dxi+dxip1
      dyp=dyj+dyjp1
      dzp=dzk+dzkp1
!
      indxp=1.D0/dxp
      indyp=1.D0/dyp
      indzp=1.D0/dzp
!
      IF ( pmodel == 1 ) THEN
         ep_e = (dxi*ep(ijke) + dxip1*ep(ijk)) * indxp
         ep_n = (dyj*ep(ijkn) + dyjp1*ep(ijk)) * indyp
         ep_t = (dzk*ep(ijkt) + dzkp1*ep(ijk)) * indzp
       ELSE IF ( pmodel == 2 ) THEN
         ep_e = 1.D0
         ep_n = 1.D0
         ep_t = 1.D0
       END IF
!
      bu(1) = rug(ijk) + alpha * dt * indxp *2.D0* ep_e * (pijk-p(ijke))
      bw(1) = rwg(ijk) + alpha * dt * indzp *2.D0* ep_t * (pijk-p(ijkt))
      bu(1) = bu(1) + dt*alphagrav*indxp*gravx*(dxi*rgp(ijke) + dxip1*rgp(ijk))
      bw(1) = bw(1) + dt*alphagrav*indzp*gravz*(dzk*rgp(ijkt) + dzkp1*rgp(ijk))
      au(1,1)=alpha*(dxi*appu(ijke,1)+dxip1*appu(ijk,1))*indxp
      au(1,1)=au(1,1)+(dxi*rgp(ijke)+dxip1*rgp(ijk))*indxp
      aw(1,1)=alpha*(dzk*appw(ijkt,1)+dzkp1*appw(ijk,1))*indzp
      aw(1,1)=aw(1,1)+(dzk*rgp(ijkt)+dzkp1*rgp(ijk))*indzp
!
      IF (job_type == JOB_TYPE_3D) THEN
        bv(1)  = rvg(ijk) + alpha * dt * indyp *2.D0* ep_n * (pijk-p(ijkn))
        bv(1) = bv(1) + dt*alphagrav*indyp*gravy*(dyj*rgp(ijkn) + dyjp1*rgp(ijk))
        av(1,1)=alpha*(dyj*appv(ijkn,1)+dyjp1*appv(ijk,1))*indyp
        av(1,1)=av(1,1)+(dyj*rgp(ijkn)+dyjp1*rgp(ijk))*indyp
      END IF
      
      DO l=2,nphase
!
! ... Explicit terms in the linear system
!
        eps_e = (dxi*rlk(ijke,l-1) + dxip1*rlk(ijk,l-1)) * indxp * inrl(l-1)
        eps_n = (dyj*rlk(ijkn,l-1) + dyjp1*rlk(ijk,l-1)) * indyp * inrl(l-1)
        eps_t = (dzk*rlk(ijkt,l-1) + dzkp1*rlk(ijk,l-1)) * indzp * inrl(l-1)
        IF ( pmodel == 1 ) THEN
           eps_e = (dxi*rlk(ijke,l-1) + dxip1*rlk(ijk,l-1)) * indxp * inrl(l-1)
           eps_n = (dyj*rlk(ijkn,l-1) + dyjp1*rlk(ijk,l-1)) * indyp * inrl(l-1)
           eps_t = (dzk*rlk(ijkt,l-1) + dzkp1*rlk(ijk,l-1)) * indzp * inrl(l-1)
        ELSE IF ( pmodel == 2 ) THEN
           eps_e = 0.D0
           eps_n = 0.D0
           eps_t = 0.D0
        END IF
        bu(l) = rus(ijk,l-1) + alpha * dt * indxp *2.D0* eps_e * (pijk-p(ijke))
        bw(l) = rws(ijk,l-1) + alpha * dt * indzp *2.D0* eps_t * (pijk-p(ijkt))

        bu(l) = bu(l) + dt*alphagrav*indxp*gravx*(dxi*rlk(ijke,l-1) + dxip1*rlk(ijk,l-1))
        bw(l) = bw(l) + dt*alphagrav*indzp*gravz*(dzk*rlk(ijkt,l-1) + dzkp1*rlk(ijk,l-1))
!
        IF (job_type == JOB_TYPE_3D) THEN
          bv(l) = rvs(ijk,l-1) + alpha * dt * indyp *2.D0* eps_n * (pijk-p(ijkn))
          bv(l) = bv(l) + dt*alphagrav*indyp*gravy*(dyj*rlk(ijkn,l-1) + dyjp1*rlk(ijk,l-1))
        END IF
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
          au(l,ll)=alpha*(dxi*appu(ijke,ls)+dxip1*appu(ijk,ls))*indxp
          au(ll,l)=au(l,ll)
          aw(l,ll)=alpha*(dzk*appw(ijkt,ls)+dzkp1*appw(ijk,ls))*indzp
          aw(ll,l)=aw(l,ll)

          IF (job_type == JOB_TYPE_3D) THEN
            av(l,ll)=alpha*(dyj*appv(ijkn,ls)+dyjp1*appv(ijk,ls))*indyp
            av(ll,l)=av(l,ll)
          END IF

        END DO
!
! ... Diagonal elements
        au(l,l)=au(l,l)+(dxi*rlk(ijke,l-1)+dxip1*rlk(ijk,l-1))*indxp
        aw(l,l)=aw(l,l)+(dzk*rlk(ijkt,l-1)+dzkp1*rlk(ijk,l-1))*indzp

        IF (job_type == JOB_TYPE_3D) av(l,l)=av(l,l)+(dyj*rlk(ijkn,l-1)+dyjp1*rlk(ijk,l-1))*indyp

      END DO
!
      RETURN
      END SUBROUTINE assemble_matrix
!----------------------------------------------------------------------
      SUBROUTINE solve_all_velocities(ijk)
! ... solve the momentum matrix an all cell faces by using
! ... Gauss-Jordan direct inversion method
!
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nphase
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE grid, ONLY: flag
      USE set_indexes, ONLY: imjk, ijmk, ijkm, ipjk, ijpk, ijkp
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ijk
!
      INTEGER :: ll, lp1, l, lj, li
      REAL*8 :: div, amul
      LOGICAL :: compute
!
! ... Use Gauss-Jordan method for matrix inversion
! ... Velocities adjacent to a boundary cell are NOT computed
!
      compute = BTEST(flag(imjk),0)
      IF (compute) THEN
        DO l=2,nphase
          IF(DABS(au1(l,l)) < plim(l-1)) THEN
            DO ll=1,nphase
              au1(ll,ll)=au1(ll,ll)+au1(l,ll)
              au1(l,ll)=0.D0
              au1(ll,l)=0.D0
            END DO
            bu1(l)=0.D0
          END IF
        END DO
        !
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
        ! ... Explicitly compute West velocities
        !
        ug(imjk) = bu1(1)
        DO l=2,nphase
          us(imjk,l-1) = bu1(l)
        END DO
        !
      END IF
!
      IF (job_type == JOB_TYPE_3D) THEN
        compute = BTEST(flag(ijmk),0)
        IF (compute) THEN
          DO l=2,nphase
            IF(DABS(av1(l,l)) < plim(l-1)) THEN
              DO ll=1,nphase
                av1(ll,ll)=av1(ll,ll)+av1(l,ll)
                av1(l,ll)=0.D0
                av1(ll,l)=0.D0
              END DO
              bv1(l)=0.D0
            END IF
          END DO
          !
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
          ! ... Explicitly compute South velocities
          !
          vg(ijmk) = bv1(1)
          DO l=2, nphase
            vs(ijmk,l-1) = bv1(l)
          END DO
          !
        END IF
      END IF
!
      compute = BTEST(flag(ijkm),0)
      IF (compute) THEN
        DO l=2,nphase
          IF(DABS(aw1(l,l)) < plim(l-1)) THEN
            DO ll=1,nphase
              aw1(ll,ll)=aw1(ll,ll)+aw1(l,ll)
              aw1(l,ll)=0.D0
              aw1(ll,l)=0.D0
            END DO
            bw1(l)=0.D0
            END IF
        END DO
        !
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
        ! ... Explicitly compute Bottom velocities
        !
        wg(ijkm) = bw1(1)
        DO l=2, nphase
          ws(ijkm,l-1) = bw1(l)
        END DO
        !
      END IF
!
      CALL solve_velocities(ijk)
!
       RETURN
       END SUBROUTINE solve_all_velocities
!----------------------------------------------------------------------
      SUBROUTINE solve_velocities(ijk)
! ... solve the momentum matrix only on East, North, and Top 
! ... cell faces by using Gauss-Jordan direct inversion method
!
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nphase
      USE domain_mapping, ONLY:  meshinds
      USE grid, ONLY: flag
      USE set_indexes, ONLY: ipjk, ijpk, ijkp
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ijk

      INTEGER :: i,j,k,imesh
      INTEGER :: ll, lp1, l, lj, li
      INTEGER :: fle, fln, flt
      REAL*8 :: div, amul
      LOGICAL :: compute
!
! ... Use Gauss-Jordan method for matrix inversion
! ... Velocities adjacent to a boundary cell are NOT computed
!
      compute = BTEST(flag(ipjk),0)
      IF (compute) THEN
        DO l=2,nphase
          IF(DABS(au(l,l)) < plim(l-1)) THEN
            DO ll=1,nphase
              au(ll,ll)=au(ll,ll)+au(l,ll)
              au(l,ll)=0.D0
              au(ll,l)=0.D0
            END DO
            bu(l)=0.D0
          END IF
        END DO
        !
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
        ! ... Explicitly compute East velocities
        !
        ug(ijk)=bu(1)
        DO l=2,nphase
          us(ijk,l-1)=bu(l)
        END DO
        !
      END IF
!
      IF (job_type == JOB_TYPE_3D) THEN
        compute = BTEST(flag(ijpk),0)
        IF (compute) THEN
          DO l=2,nphase
            IF(DABS(av(l,l)) < plim(l-1)) THEN
              DO ll=1,nphase
                av(ll,ll)=av(ll,ll)+av(l,ll)
                av(l,ll)=0.D0
                av(ll,l)=0.D0
              END DO
              bv(l)=0.D0
            END IF
          END DO
          !
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
          ! ... Explicitly compute North velocities
          !
          vg(ijk)=bv(1)
          DO l=2,nphase
            vs(ijk,l-1)=bv(l)
          END DO
          !
        END IF
      END IF
!
      compute = BTEST(flag(ijkp),0)
      IF (compute) THEN
        DO l=2,nphase
          IF(DABS(aw(l,l)) < plim(l-1)) THEN
            DO ll=1,nphase
              aw(ll,ll)=aw(ll,ll)+aw(l,ll)
              aw(l,ll)=0.D0
              aw(ll,l)=0.D0
            END DO
            bw(l)=0.D0
          END IF
        END DO
        !
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
        ! ... Explicitly compute Top velocities
        !
        wg(ijk)=bw(1)
        DO l=2,nphase
          ws(ijk,l-1)=bw(l)
        END DO
        !
      END IF
!
      RETURN
      END SUBROUTINE solve_velocities
!----------------------------------------------------------------------
!     H E R E   S T A R T   T H E   O P T I M I Z E D   R O U T I N E S 
!----------------------------------------------------------------------
      SUBROUTINE matspre_3phase( au1t, aut, av1t, avt, aw1t, awt, i, j, k, ijk )
!
! ... Computes matrix elements to solve momentum-balance 
! ... linear system of coupled equations in every cell face
! ... OPTIMIZED ROUTINE FOR THREE-DIMENSIONAL THREE PHASE SIMULATIONS
!
      USE grid, ONLY: dx, dy, dz
      USE set_indexes, ONLY: ijke, ijkn, ijkt, ijkw, ijks, ijkb
      USE tilde_momentum, ONLY: appu, appv,  appw

      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ijk, i, j, k
      REAL*8 :: au1t(6), aut(6)
      REAL*8 :: av1t(6), avt(6)
      REAL*8 :: aw1t(6), awt(6)
!
      REAL*8 :: indxp, indyp, indzp, indxm, indym, indzm
      REAL*8 :: dxdm, dxmdm, dxdp, dxpdp
      REAL*8 :: dydm, dymdm, dydp, dypdp
      REAL*8 :: dzdm, dzmdm, dzdp, dzpdp
!
      indxm=1.D0/( dx(i)+dx(i-1) )
      indxp=1.D0/( dx(i)+dx(i+1) )
      dxdm = dx(i) * indxm
      dxmdm= dx(i-1) * indxm
      dxdp = dx(i) * indxp
      dxpdp= dx(i+1) * indxp

      indym=1.D0/( dy(j)+dy(j-1) )
      indyp=1.D0/( dy(j)+dy(j+1) )
      dydm = dy(j) * indym
      dymdm= dy(j-1) * indym
      dydp = dy(j) * indyp
      dypdp= dy(j+1) * indyp

      indzm=1.D0/( dz(k)+dz(k-1) )
      indzp=1.D0/( dz(k)+dz(k+1) )
      dzdm = dz(k) * indzm
      dzmdm= dz(k-1) * indzm
      dzdp = dz(k) * indzp
      dzpdp= dz(k+1) * indzp

      au1t(1) = ( dxdm * appu(ijkw,1) + dxmdm * appu(ijk,1) )
      aut(1)  = ( dxdp * appu(ijke,1) + dxpdp * appu(ijk,1) )
      au1t(2) = ( dxdm * appu(ijkw,2) + dxmdm * appu(ijk,2) )
      aut(2)  = ( dxdp * appu(ijke,2) + dxpdp * appu(ijk,2) )
      au1t(3) = ( dxdm * appu(ijkw,3) + dxmdm * appu(ijk,3) )
      aut(3)  = ( dxdp * appu(ijke,3) + dxpdp * appu(ijk,3) )
      au1t(4) = ( dxdm * appu(ijkw,4) + dxmdm * appu(ijk,4) )
      aut(4)  = ( dxdp * appu(ijke,4) + dxpdp * appu(ijk,4) )
      au1t(5) = ( dxdm * appu(ijkw,5) + dxmdm * appu(ijk,5) )
      aut(5)  = ( dxdp * appu(ijke,5) + dxpdp * appu(ijk,5) )
      au1t(6) = ( dxdm * appu(ijkw,6) + dxmdm * appu(ijk,6) )
      aut(6)  = ( dxdp * appu(ijke,6) + dxpdp * appu(ijk,6) )

      av1t(1) = ( dydm * appv(ijks,1) + dymdm * appv(ijk,1) )
      avt(1)  = ( dydp * appv(ijkn,1) + dypdp * appv(ijk,1) )
      av1t(2) = ( dydm * appv(ijks,2) + dymdm * appv(ijk,2) )
      avt(2)  = ( dydp * appv(ijkn,2) + dypdp * appv(ijk,2) )
      av1t(3) = ( dydm * appv(ijks,3) + dymdm * appv(ijk,3) )
      avt(3)  = ( dydp * appv(ijkn,3) + dypdp * appv(ijk,3) )
      av1t(4) = ( dydm * appv(ijks,4) + dymdm * appv(ijk,4) )
      avt(4)  = ( dydp * appv(ijkn,4) + dypdp * appv(ijk,4) )
      av1t(5) = ( dydm * appv(ijks,5) + dymdm * appv(ijk,5) )
      avt(5)  = ( dydp * appv(ijkn,5) + dypdp * appv(ijk,5) )
      av1t(6) = ( dydm * appv(ijks,6) + dymdm * appv(ijk,6) )
      avt(6)  = ( dydp * appv(ijkn,6) + dypdp * appv(ijk,6) )

      aw1t(1) = ( dzdm * appw(ijkb,1) + dzmdm * appw(ijk,1) )
      awt(1)  = ( dzdp * appw(ijkt,1) + dzpdp * appw(ijk,1) )
      aw1t(2) = ( dzdm * appw(ijkb,2) + dzmdm * appw(ijk,2) )
      awt(2)  = ( dzdp * appw(ijkt,2) + dzpdp * appw(ijk,2) )
      aw1t(3) = ( dzdm * appw(ijkb,3) + dzmdm * appw(ijk,3) )
      awt(3)  = ( dzdp * appw(ijkt,3) + dzpdp * appw(ijk,3) )
      aw1t(4) = ( dzdm * appw(ijkb,4) + dzmdm * appw(ijk,4) )
      awt(4)  = ( dzdp * appw(ijkt,4) + dzpdp * appw(ijk,4) )
      aw1t(5) = ( dzdm * appw(ijkb,5) + dzmdm * appw(ijk,5) )
      awt(5)  = ( dzdp * appw(ijkt,5) + dzpdp * appw(ijk,5) )
      aw1t(6) = ( dzdm * appw(ijkb,6) + dzmdm * appw(ijk,6) )
      awt(6)  = ( dzdp * appw(ijkt,6) + dzpdp * appw(ijk,6) )

      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------
      SUBROUTINE matsvels_3phase( au1t, aut, av1t, avt, aw1t, awt, i, j, k, ijk )
!
! ... Computes matrix elements to solve momentum-balance 
! ... linear system of coupled equations in every cell face
! ... OPTIMIZED ROUTINE FOR THREE-DIMENSIONAL THREE PHASE SIMULATIONS
!
      USE gas_solid_density, ONLY: rgp, rlk
      USE grid, ONLY: dx, dy, dz
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p, ep, pmodel
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes, ONLY: ipjk, ijpk, ijkp
      USE set_indexes, ONLY: ijke, ijkn, ijkt, ijkw, ijks, ijkb
      USE tilde_momentum, ONLY: rug, rvg, rwg, rus, rvs, rws
      USE time_parameters, ONLY: dt
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE grid, ONLY: flag


      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ijk, i, j, k
      REAL*8 :: au1t(6), aut(6)
      REAL*8 :: av1t(6), avt(6)
      REAL*8 :: aw1t(6), awt(6)
!
      INTEGER, SAVE :: ijk_old = -1
      REAL*8, SAVE :: dxp, dyp, dzp, dxm, dym, dzm
      REAL*8, SAVE :: indxp, indyp, indzp, indxm, indym, indzm
      REAL*8, SAVE :: dxdm, dxmdm, dxdp, dxpdp
      REAL*8, SAVE :: dydm, dymdm, dydp, dypdp
      REAL*8, SAVE :: dzdm, dzmdm, dzdp, dzpdp
      REAL*8, SAVE :: pw, pe, ps, pn, pt, pb
      REAL*8, SAVE :: depw, depe, deps, depn, dept, depb
      REAL*8, SAVE :: drgpw, drgpe, drgps, drgpn, drgpt, drgpb
      REAL*8, SAVE :: drlkw1, drlkw2
      REAL*8, SAVE :: drlks1, drlks2
      REAL*8, SAVE :: drlkb1, drlkb2
      REAL*8, SAVE :: drlke1, drlke2
      REAL*8, SAVE :: drlkn1, drlkn2
      REAL*8, SAVE :: drlkt1, drlkt2
      REAL*8, SAVE :: epijk, rgpijk, rlk1, rlk2
      REAL*8, SAVE :: pijk, dpijk
      REAL*8, SAVE :: twodt
      REAL*8, SAVE :: div, amul
      LOGICAL, SAVE :: flim, flip, fljm, fljp, flkm, flkp
      REAL*8 :: au1(3,3), bu1(3), av1(3,3), bv1(3), aw1(3,3), bw1(3)
      REAL*8 :: au(3,3), bu(3), av(3,3), bv(3), aw(3,3), bw(3)
!
      IF( ijk /= ijk_old ) THEN

        dxm  =dx(i)+dx(i-1)
        dxp  =dx(i)+dx(i+1)
        indxm=1.D0/dxm
        indxp=1.D0/dxp
        dxdm =dx(i)   * indxm
        dxmdm=dx(i-1) * indxm
        dxdp =dx(i)   * indxp
        dxpdp=dx(i+1) * indxp

        dym  =dy(j)+dy(j-1)
        dyp  =dy(j)+dy(j+1)
        indym=1.D0/dym
        indyp=1.D0/dyp
        dydm =dy(j)   * indym
        dymdm=dy(j-1) * indym
        dydp =dy(j)   * indyp
        dypdp=dy(j+1) * indyp

        dzm  =dz(k)+dz(k-1)
        dzp  =dz(k)+dz(k+1)
        indzm=1.D0/dzm
        indzp=1.D0/dzp
        dzdm =dz(k)   * indzm
        dzmdm=dz(k-1) * indzm
        dzdp =dz(k)   * indzp
        dzpdp=dz(k+1) * indzp

        flim = BTEST( flag(imjk), 0 )
        fljm = BTEST( flag(ijmk), 0 )
        flkm = BTEST( flag(ijkm), 0 )
        flip = BTEST( flag(ipjk), 0 )
        fljp = BTEST( flag(ijpk), 0 )
        flkp = BTEST( flag(ijkp), 0 )

        pw   = p(ijkw)
        pe   = p(ijke)
        ps   = p(ijks)
        pn   = p(ijkn)
        pb   = p(ijkb)
        pt   = p(ijkt)

        depw   = dxdm * ep(ijkw)
        depe   = dxdp * ep(ijke)
        deps   = dydm * ep(ijks)
        depn   = dydp * ep(ijkn)
        depb   = dzdm * ep(ijkb)
        dept   = dzdp * ep(ijkt)

        drgpw  = dxdm * rgp(ijkw)
        drgpe  = dxdp * rgp(ijke)
        drgps  = dydm * rgp(ijks)
        drgpn  = dydp * rgp(ijkn)
        drgpb  = dzdm * rgp(ijkb)
        drgpt  = dzdp * rgp(ijkt)

        drlkw1 = dxdm * rlk(ijkw,1) 
        drlks1 = dydm * rlk(ijks,1) 
        drlkb1 = dzdm * rlk(ijkb,1) 
        drlke1 = dxdp * rlk(ijke,1) 
        drlkn1 = dydp * rlk(ijkn,1) 
        drlkt1 = dzdp * rlk(ijkt,1) 

        drlkw2 = dxdm * rlk(ijkw,2) 
        drlks2 = dydm * rlk(ijks,2) 
        drlkb2 = dzdm * rlk(ijkb,2) 
        drlke2 = dxdp * rlk(ijke,2) 
        drlkn2 = dydp * rlk(ijkn,2) 
        drlkt2 = dzdp * rlk(ijkt,2) 

        ijk_old = ijk

      END IF

      pijk    = p(ijk)
      epijk   = ep(ijk)
      rgpijk  = rgp(ijk)
      rlk2    = rlk(ijk,2)
      rlk1    = rlk(ijk,1)
!
      twodt = 2.0d0 * dt

!
! ... solve the momentum matrix an all cell faces by using
! ... Gauss-Jordan direct inversion method
!
!
      ! ... Use Gauss-Jordan method for matrix inversion

      au1(1,1) = au1t(1) + ( drgpw + dxmdm * rgpijk )
      au1(2,1) = au1t(2)
      au1(3,1) = au1t(4)
      au1(1,2) = au1t(2)
      au1(2,2) = au1t(3) + ( drlkw1 + dxmdm * rlk1 )
      au1(3,2) = au1t(5)
      au1(1,3) = au1t(4)
      au1(2,3) = au1t(5)
      au1(3,3) = au1t(6) + ( drlkw2 + dxmdm * rlk2 )

      bu1(1) = ( depw + dxmdm * epijk )
      bu1(2) = ( drlkw1 + dxmdm * rlk1 ) * inrl(1)
      bu1(3) = ( drlkw2 + dxmdm * rlk2 ) * inrl(2)
      dpijk  = ( pw - pijk ) * indxm * twodt
      bu1(1) = rug(imjk)   + bu1(1) * dpijk
      bu1(2) = rus(imjk,1) + bu1(2) * dpijk
      bu1(3) = rus(imjk,2) + bu1(3) * dpijk
!
      IF ( flim ) THEN

        call matsvels_solve( au1, bu1 )

        ug(imjk)   = bu1(1)
        us(imjk,1) = bu1(2)
        us(imjk,2) = bu1(3)

      END IF


      av1(1,1) = av1t(1) + ( drgps + dymdm * rgpijk )
      av1(2,1) = av1t(2)
      av1(3,1) = av1t(4)
      av1(1,2) = av1t(2)
      av1(2,2) = av1t(3) + ( drlks1 + dymdm * rlk1 )
      av1(3,2) = av1t(5)
      av1(1,3) = av1t(4)
      av1(2,3) = av1t(5)
      av1(3,3) = av1t(6) + ( drlks2 + dymdm * rlk2 )

      bv1(1) = ( deps + dymdm * epijk )
      bv1(2) = ( drlks1 + dymdm * rlk1 ) * inrl(1)
      bv1(3) = ( drlks2 + dymdm * rlk2 ) * inrl(2)
      dpijk  = ( ps - pijk ) * indym * twodt
      bv1(1) = rvg(ijmk)   + bv1(1) * dpijk
      bv1(2) = rvs(ijmk,1) + bv1(2) * dpijk
      bv1(3) = rvs(ijmk,2) + bv1(3) * dpijk

      IF ( fljm ) THEN

          call matsvels_solve( av1, bv1 )
!
          vg(ijmk)   = bv1(1)
          vs(ijmk,1) = bv1(2)
          vs(ijmk,2) = bv1(3)
!
      END IF

      aw1(1,1) = aw1t(1) + ( drgpb + dzmdm * rgpijk )
      aw1(2,1) = aw1t(2)
      aw1(3,1) = aw1t(4)
      aw1(1,2) = aw1t(2)
      aw1(2,2) = aw1t(3) + ( drlkb1 + dzmdm * rlk1 )
      aw1(3,2) = aw1t(5)
      aw1(1,3) = aw1t(4)
      aw1(2,3) = aw1t(5)
      aw1(3,3) = aw1t(6) + ( drlkb2 + dzmdm * rlk2 )

      bw1(1) = ( depb + dzmdm * epijk )
      bw1(2) = ( drlkb1 + dzmdm * rlk1 ) * inrl(1)
      bw1(3) = ( drlkb2 + dzmdm * rlk2 ) * inrl(2)
      dpijk  = ( pb - pijk ) * indzm * twodt
      bw1(1) = rwg(ijkm)   + bw1(1) * dpijk
      bw1(2) = rws(ijkm,1) + bw1(2) * dpijk
      bw1(3) = rws(ijkm,2) + bw1(3) * dpijk


      IF ( flkm ) THEN

        call matsvels_solve( aw1, bw1 )
!
        wg(ijkm) = bw1(1)
        ws(ijkm,2-1) = bw1(2)
        ws(ijkm,3-1) = bw1(3)
!
      END IF
!
      ! ... solve the momentum matrix only on East, North, and Top 
      ! ... cell faces by using Gauss-Jordan direct inversion method
!
      ! ... Use Gauss-Jordan method for matrix inversion

      au(1,1)  = aut(1)  + ( drgpe + dxpdp * rgpijk )
      au(2,1)  = aut(2)
      au(3,1)  = aut(4)
      au(1,2)  = aut(2)
      au(2,2)  = aut(3)  + ( drlke1 + dxpdp * rlk1 )
      au(3,2)  = aut(5)
      au(1,3)  = aut(4)
      au(2,3)  = aut(5)
      au(3,3)  = aut(6)  + ( drlke2 + dxpdp * rlk2 )

      bu(1) = ( depe + dxpdp * epijk )
      bu(2) = ( drlke1 + dxpdp * rlk1 ) * inrl(1)
      bu(3) = ( drlke2 + dxpdp * rlk2 ) * inrl(2)
      dpijk  = ( pijk - pe ) * indxp * twodt
      bu(1)  = rug(ijk)    + bu(1) * dpijk
      bu(2)  = rus(ijk,1)  + bu(2) * dpijk
      bu(3)  = rus(ijk,2)  + bu(3) * dpijk
!

      IF ( flip ) THEN

        call matsvels_solve( au, bu )
!
        ug(ijk)=bu(1)
        us(ijk,1)=bu(2)
        us(ijk,2)=bu(3)
!
      END IF


      av(1,1)  = avt(1)  + ( drgpn + dypdp * rgpijk )
      av(2,1)  = avt(2)
      av(3,1)  = avt(4)
      av(1,2)  = avt(2)
      av(2,2)  = avt(3)  + ( drlkn1 + dypdp * rlk1 )
      av(3,2)  = avt(5)
      av(1,3)  = avt(4)
      av(2,3)  = avt(5)
      av(3,3)  = avt(6)  + ( drlkn2 + dypdp * rlk2 )

      bv(1) = ( depn + dypdp * epijk )
      bv(2) = ( drlkn1 + dypdp * rlk1 ) * inrl(1)
      bv(3) = ( drlkn2 + dypdp * rlk2 ) * inrl(2)
      dpijk  = ( pijk - pn ) * indyp * twodt
      bv(1)  = rvg(ijk)    + bv(1) * dpijk
      bv(2)  = rvs(ijk,1)  + bv(2) * dpijk
      bv(3)  = rvs(ijk,2)  + bv(3) * dpijk

      IF ( fljp ) THEN

         call matsvels_solve( av, bv )
!
         vg(ijk)  = bv(1)
         vs(ijk,1)= bv(2)
         vs(ijk,2)= bv(3)
!
      END IF
!
      aw(1,1)  = awt(1)  + ( drgpt + dzpdp * rgpijk )
      aw(2,1)  = awt(2)
      aw(3,1)  = awt(4)
      aw(1,2)  = awt(2)
      aw(2,2)  = awt(3)  + ( drlkt1 + dzpdp * rlk1 )
      aw(3,2)  = awt(5)
      aw(1,3)  = awt(4)
      aw(2,3)  = awt(5)
      aw(3,3)  = awt(6)  + ( drlkt2 + dzpdp * rlk2 )

      bw(1) = ( dept   + dzpdp * epijk )
      bw(2) = ( drlkt1 + dzpdp * rlk1 ) * inrl(1)
      bw(3) = ( drlkt2 + dzpdp * rlk2 ) * inrl(2)
      dpijk  = ( pijk - pt ) * indzp * twodt
      bw(1)  = rwg(ijk)    + bw(1) * dpijk
      bw(2)  = rws(ijk,1)  + bw(2) * dpijk
      bw(3)  = rws(ijk,2)  + bw(3) * dpijk

      IF ( flkp ) THEN

        call matsvels_solve( aw, bw )
!
        wg(ijk)   = bw(1)
        ws(ijk,1) = bw(2)
        ws(ijk,2) = bw(3)
!
      END IF
!
      RETURN

      END SUBROUTINE matsvels_3phase
!----------------------------------------------------------------------
      SUBROUTINE matsvels_solve( a, b )

      REAL*8 :: a(3,3), b(3)
      REAL*8 :: d, x1, x2, amul

      IF( DABS(a(2,2)) < plim(1) ) THEN

        a(1,1) = a(1,1) + a(2,1)
        a(3,3) = a(3,3) + a(3,2)

        IF( DABS(a(3,3)) < plim(2) ) THEN

          a(1,1) = a(1,1) + a(3,1)
          a(2,2) = a(2,2) + a(3,2)

          IF( a(1,1) /= 0 ) THEN

            !write(6,*) 'c111 ok'

            b(1) = b(1) / a(1,1)
            b(2) = 0.D0
            b(3) = 0.D0

          ELSE

            !write(6,*) 'c112 ok'

            b(1) = b(1)
            b(2) = 0.D0
            b(3) = 0.D0

          END IF

        ELSE

          IF( a(1,1) /= 0 ) THEN

            !write(6,*) 'c121 ok'  

            d      = 1.D0 / a(1,1)
            a(1,3) = a(1,3) * d
            b(1)   = b(1)   * d
            ! a(1,1) = 0.D0
            amul   = a(3,1)
            a(3,3) = a(3,3) - amul * a(1,3)
            b(3)   = b(3)   - amul * b(1)

            d      = 1.D0 / a(3,3)
            b(3)   = b(3) * d
            ! a(3,3) = 0.D0
            amul   = a(1,3)
            b(1)   = b(1) - amul * b(3)
    
            b(2)   = 0.0d0

            !d = ( a(1,1) * a(3,3) - a(1,3) * a(3,1) )
            !d = 1.0d0 / d
            !x1 = (  a(3,3) * b(1) - a(3,1) * b(3) ) * d
            !x2 = ( -a(1,3) * b(1) + a(1,1) * b(3) ) * d
            !b(1) = x1
            !b(2) = 0.0d0
            !b(3) = x2
        
          ELSE
    
            !write(6,*) 'c122 ok'
    
            d      = 1.D0 / a(3,3)
            b(3)   = b(3) * d
            ! a(3,3) = 0.D0
            amul   = a(1,3)
            b(1)   = b(1) - amul * b(3)
    
            b(2)   = 0.0d0
    
          END IF
    
        END IF
    
      ELSE IF ( DABS(a(3,3)) < plim(2) ) THEN
    
        a(1,1) = a(1,1) + a(3,1)
        a(2,2) = a(2,2) + a(3,2)

        IF( a(1,1) /= 0 ) THEN

            !write(6,*) 'c21 ok'

            d      = 1.D0 / a(1,1)
            a(1,2) = a(1,2) * d
            b(1)   = b(1)   * d
            ! a(1,1) = 0.D0
            amul   = a(2,1)
            a(2,2) = a(2,2) - amul * a(1,2)
            b(2)   = b(2)   - amul * b(1)
    
            d      = 1.D0/a(2,2)
            b(2)   = b(2)*d
            ! a(2,2) = 0.D0
            amul   = a(1,2)
            b(1)   = b(1)-amul*b(2)
    
    
            !d = ( a(1,1) * a(2,2) - a(1,2) * a(2,1) )
            !d = 1.0d0 / d
            !x1 = (  a(2,2) * b(1) - a(2,1) * b(2) ) * d
            !x2 = ( -a(1,2) * b(1) + a(1,1) * b(2) ) * d
            !b(1) = x1
            !b(2) = x2
    
            b(3) = 0.0d0
    
        ELSE
    
            !write(6,*) 'c22 ok'

            d      = 1.D0 / a(2,2)
            b(2)   = b(2) * d
            ! a(2,2) = 0.D0
            amul   = a(1,2)
            b(1)   = b(1) - amul * b(2)
    
            b(3) = 0.0d0
    
        END IF

      ELSE

        IF( a(1,1) /= 0.D0 ) THEN

            !write(6,*) 'c31'

            d      = 1.D0 / a(1,1)
            a(1,2) = a(1,2) * d
            a(1,3) = a(1,3) * d
            b(1)   = b(1)   * d
            ! a(1,1) = 0.D0
            amul   = a(2,1)
            a(2,2) = a(2,2) - amul * a(1,2)
            a(2,3) = a(2,3) - amul * a(1,3)
            b(2)   = b(2)   - amul * b(1)
            amul   = a(3,1)
            a(3,2) = a(3,2) - amul * a(1,2)
            a(3,3) = a(3,3) - amul * a(1,3)
            b(3)   = b(3)   - amul * b(1)

            d      = 1.D0 / a(2,2)
            a(2,3) = a(2,3) * d
            b(2)   = b(2)   * d
            ! a(2,2) = 0.D0
            amul   = a(1,2)
            a(1,3) = a(1,3) - amul * a(2,3)
            b(1)   = b(1)   - amul * b(2)
            amul   = a(3,2)
            a(3,3) = a(3,3) - amul * a(2,3)
            b(3)   = b(3)   - amul * b(2)

            d      = 1.D0 / a(3,3)
            b(3)   = b(3) * d
            ! a(3,3) = 0.D0
            amul   = a(1,3)
            b(1)   = b(1) - amul * b(3)
            amul   = a(2,3)
            b(2)   = b(2) - amul * b(3)

        ELSE

            !write(6,*) 'c32'

            !d = ( a(2,2) * a(3,3) - a(2,3) * a(3,2) )
            !d = 1.0d0 / d
            !x1 = (  a(3,3) * b(2) - a(3,2) * b(3) ) * d
            !x2 = ( -a(2,3) * b(2) + a(2,2) * b(3) ) * d
            !b(1) = 0.0d0
            !b(2) = x1
            !b(3) = x2

            d      = 1.D0 / a(2,2)
            a(2,3) = a(2,3) * d
            b(2)   = b(2)   * d
            ! a(2,2) = 0.D0
            amul   = a(1,2)
            a(1,3) = a(1,3) - amul * a(2,3)
            b(1)   = b(1)   - amul * b(2)
            amul   = a(3,2)
            a(3,3) = a(3,3) - amul * a(2,3)
            b(3)   = b(3)   - amul * b(2)
    
            d      = 1.D0 / a(3,3)
            b(3)   = b(3) * d
            ! a(3,3) = 0.D0
            amul   = a(1,3)
            b(1)   = b(1) - amul * b(3)
            amul   = a(2,3)
            b(2)   = b(2) - amul * b(3)
    
        END IF

      END IF

      END SUBROUTINE matsvels_solve
!----------------------------------------------------------------------
      SUBROUTINE matsvels_solve_old( a, b )

      real*8 :: a(3,3), b(3)
      real*8 :: div, amul

      IF( DABS(a(2,2)) < plim(1) ) THEN

              !write(6,*) 'if1'

              a(1,1) = a(1,1) + a(2,1)
              a(2,1) = 0.D0
              a(1,2) = 0.D0

              a(2,2) = 0.D0

              a(3,3) = a(3,3) + a(2,3)
              a(2,3) = 0.D0
              a(3,2) = 0.D0

              b(2) = 0.D0

      END IF


      IF( DABS(a(3,3)) < plim(2) ) THEN

              !write(6,*) 'if2'

              a(1,1) = a(1,1) + a(3,1)
              a(3,1) = 0.D0
              a(1,3) = 0.D0

              a(2,2) = a(2,2) + a(3,2)
              a(3,2) = 0.D0
              a(2,3) = 0.D0

              a(3,3) = 0.D0

              b(3) = 0.D0

      END IF

      IF( a(1,1) /= 0.D0 ) THEN

            !write(6,*) 'if3'

            div = 1.D0 / a(1,1)
            a(1,2) = a(1,2) * div
            a(1,3) = a(1,3) * div
            b(1)   = b(1)   * div
            a(1,1) = 0.D0
            !li=2
              amul = a(2,1)
              a(2,2) = a(2,2) - amul * a(1,2)
              a(2,3) = a(2,3) - amul * a(1,3)
              b(2)   = b(2)   - amul * b(1)
            !li=3
              amul = a(3,1)
              a(3,2) = a(3,2) - amul * a(1,2)
              a(3,3) = a(3,3) - amul * a(1,3)
              b(3)   = b(3)   - amul * b(1)

      END IF

      IF( a(2,2) /= 0.D0 ) THEN
            !write(6,*) 'if4'
            div=1.D0/a(2,2)
            a(2,3)=a(2,3)*div
            b(2)=b(2)*div
            a(2,2)=0.D0
            !li=1
              amul=a(1,2)
              a(1,3)=a(1,3)-amul*a(2,3)
              b(1)=b(1)-amul*b(2)
            !li=3
              amul=a(3,2)
              a(3,3)=a(3,3)-amul*a(2,3)
              b(3)=b(3)-amul*b(2)
      END IF

      IF( a(3,3) /= 0.D0 ) THEN
            !write(6,*) 'if5'
            div=1.D0/a(3,3)
            b(3)=b(3)*div
            a(3,3)=0.D0
            !li=1
              amul=a(1,3)
              b(1)=b(1)-amul*b(3)
            !li=2
              amul=a(2,3)
              b(2)=b(2)-amul*b(3)
      END IF

      END SUBROUTINE matsvels_solve_old

!------------------------------------------------------------------------
      END MODULE phases_matrix
!------------------------------------------------------------------------
