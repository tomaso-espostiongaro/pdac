!----------------------------------------------------------------------
      MODULE enthalpy_matrix
!----------------------------------------------------------------------
      IMPLICIT NONE

! ... Interphase enthalpy matrix elements
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: bt
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: at
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE ftem
!----------------------------------------------------------------------
! ... This routine computes the matrix elements for the thermal 
! ... interphase coupling, and solves for the enthalpies
!
      USE control_flags, ONLY: job_type
      USE dimensions
      USE domain_decomposition, ONLY: ncint, meshinds
      USE eos_gas, ONLY: cg
      USE gas_solid_density, ONLY: rgp, rgpn, rlk, rlkn, rog
      USE gas_solid_temperature, ONLY: sieg, siegn, sies, siesn
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_viscosity, ONLY: mug, kapg
      USE grid, ONLY: fl_l
      USE grid, ONLY: dx, dy, dz
      USE specific_heat_module, ONLY: ck
      USE heat_transfer, ONLY: hvs
      USE pressure_epsilon, ONLY: ep, p, pn
      USE reactions, ONLY: hrex, irex
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm
      USE set_indexes, ONLY: ijke, ijkn, ijkt, ijkw, ijks, ijkb
      USE tilde_energy, ONLY: rhg, rhs
      USE time_parameters, ONLY: time, dt
      IMPLICIT NONE
!
      REAL*8 :: hv
      REAL*8 :: hrexs, hrexg
      REAL*8 :: dugs, dvgs, dwgs
      REAL*8 :: dpxyz, deltap
      REAL*8 :: indxc, indyc, indzc
      REAL*8 :: ugc, vgc, wgc
      INTEGER :: is, is1
      INTEGER :: ijk, i, j, k, imesh
!
      ALLOCATE(at(nphase, nphase))
      ALLOCATE(bt(nphase))
!
          at = 0.D0
          bt = 0.D0
          hv = 0.D0
!
!pdac------------
! Control heat reactions
      hrexg=0.D0
      hrexs=0.D0
!pdac------------
!
      indxc = 0.D0 ; indyc = 0.D0 ; indzc = 0.D0
      ugc = 0.D0 ; vgc = 0.D0 ; wgc = 0.D0
!
      DO ijk = 1, ncint

        IF(fl_l(ijk) == 1) THEN
          CALL meshinds(ijk,imesh,i,j,k)
          CALL subscr(ijk)
!
          IF(irex == 2) CALL hrex(ijk,hrexg,hrexs)
!
! ... The pressure contribution can be treated implicitly
!
          indxc = 1.D0/(dx(i)+(dx(i+1)+dx(i-1))*0.5D0)
          indzc = 1.D0/(dz(k)+(dz(k+1)+dz(k-1))*0.5D0)
          ugc = (ug(ijk)+ug(imjk))/2.D0
          wgc = (wg(ijk)+wg(ijkm))/2.D0

          IF (job_type == '3D') THEN
            indyc = 1.D0 / (dy(j)+(dy(j+1)+dy(j-1))*0.5D0)
            vgc = (vg(ijk)+vg(ijmk))/2.D0
          END IF
!
          dpxyz= dt * indxc * ugc * (p(ijke)-p(ijkw)) +   &
                 dt * indyc * vgc * (p(ijkn)-p(ijks)) +   &
                 dt * indzc * wgc * (p(ijkt)-p(ijkb))
!
          deltap = ep(ijk) * (p(ijk) - pn(ijk) + dpxyz)
!
          at(1,1) = rgp(ijk)
          bt(1)   = siegn(ijk) * rgpn(ijk) + rhg(ijk) + deltap - hrexg

          DO is=1, nsolid
            is1=is+1
!
! ... Compute gas-particle heat transfer coefficients
!
            dugs = ( (ug(ijk)-us(ijk,is)) + (ug(imjk)-us(imjk,is)) ) * 0.5D0
            dwgs = ( (wg(ijk)-ws(ijk,is)) + (wg(ijkm)-ws(ijkm,is)) ) * 0.5D0

            IF (job_type == '2D') THEN
              dvgs = 0.D0
            ELSE IF (job_type == '3D') THEN
              dvgs = ( (vg(ijk)-vs(ijk,is)) + (vg(ijmk)-vs(ijmk,is)) ) * 0.5D0
            END IF

            CALL hvs(hv, rlk(ijk,is), rog(ijk), ep(ijk), &
                     dugs, dvgs, dwgs, mug(ijk), kapg(ijk), cg(ijk), is)
!
            at(1,1)     = at(1,1)     + dt * hv / cg(ijk)
            at(1,is1)   =             - dt * hv / ck(is,ijk)
            at(is1,1)   =             - dt * hv / cg(ijk)
            at(is1,is1) = rlk(ijk,is) + dt * hv / ck(is,ijk)
!
            bt(is1) = rlkn(ijk,is) * siesn(ijk,is) + rhs(ijk, is)
          END DO
!
! ... Solve the interphase enthalpy matrix by using Gauss inversion
!
          CALL invdm(at, bt, ijk)
!
          sieg(ijk) = bt(1)
          DO is=1, nsolid
            sies(ijk,is) = bt(is+1)
          END DO
!
        END IF
      END DO
!
      DEALLOCATE(at)
      DEALLOCATE(bt)
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE invdm(a, b, ijk)
!----------------------------------------------------------------------
! ... this routine solves the N-Phases Energy-Equation Matrix
! ... (Gauss elimination)
!
      USE dimensions
      USE gas_solid_density, ONLY: rlk
      USE particles_constants, ONLY: inrl
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ijk
      REAL*8 :: a(nphase,nphase),b(nphase)
!
      INTEGER :: is 
      REAL*8 :: div
!
      DO is=nphase,2,-1
!        IF(abs(a(is,is)) <= 1.D-6) THEN
        IF(rlk(ijk,is-1)*inrl(is-1) <= 1.D-9) THEN
          a(1,is)=0.D0
          a(is,1)=0.D0
          b(is)=0.D0
        ELSE
!
! ... eliminate all cross elements in the gas enthalpy equation
!
          div=1.D0/a(is,is)
          a(is,1)=a(is,1)*div
          b(is)=b(is)*div
          b(1)=b(1)-a(1,is)*b(is)
          a(1,1)=a(1,1)-a(1,is)*a(is,1)
        ENDIF
      END DO
!
      b(1)=b(1)/a(1,1)
      DO is=2,nphase
        b(is)=b(is)-a(is,1)*b(1)
      END DO
!
      RETURN
      END SUBROUTINE invdm
!----------------------------------------------------------------------
      END MODULE enthalpy_matrix
!----------------------------------------------------------------------
