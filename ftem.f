!----------------------------------------------------------------------
      MODULE enthalpy_matrix
!----------------------------------------------------------------------
      IMPLICIT NONE

! ... Interphase enthalpy matrix elements
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: bt
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: at

      REAL*8 :: flim
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
      USE eos_gas, ONLY: cg, ygc, caloric_eosg, thermal_eosg
      USE eos_solid, ONLY: caloric_eosl
      USE gas_solid_density, ONLY: rgp, rgpn, rlk, rlkn, rog
      USE gas_solid_temperature, ONLY: sieg, siegn, sies, siesn
      USE gas_solid_temperature, ONLY: tg, ts
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_viscosity, ONLY: mug, kapg
      USE grid, ONLY: flag
      USE grid, ONLY: dx, dy, dz
      USE immersed_boundaries, ONLY: numx, numy, numz, immb, faces
      USE specific_heat_module, ONLY: ck, cp
      USE heat_transfer, ONLY: hvs
      USE particles_constants, ONLY: cps
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
      INTEGER :: ijk, i, j, k, imesh, info, num
      INTEGER :: b_e, b_w, b_t, b_b, b_n, b_s
      INTEGER :: fx, fy, fz
      REAL*8 :: ivf

      LOGICAL :: forced = .FALSE.
      info = 0
      num = 0
      fx = 0
      fy = 0
      fz = 0
!
! ... Initialize the cell fractions for immersed boundaries
!
      b_e = 1; b_w = 1; b_t = 1; b_b = 1; b_n = 1; b_s = 1; ivf = 1.D0
!
! ... Allocate and initialize the enthalpy matrix elements
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

        IF (immb == 1) THEN
          fx = numx(ijk)
          IF (job_type == '3D') fy = numy(ijk)
          fz = numz(ijk)
          forced = (fx/=0 .OR. fy/=0 .OR. fz/=0)
          CALL faces(ijk, b_e, b_w, b_t, b_b, b_n, b_s, ivf)
        END IF

        IF (flag(ijk) == 1) THEN

          CALL meshinds(ijk,imesh,i,j,k)
          CALL subscr(ijk)
!
          IF(irex == 2) CALL hrex(ijk,hrexg,hrexs)
!
! ... The pressure terms can be treated implicitly
! ... (coupled with momentum equation into the iterative solver)
! ... This term is not modified (?) by the forcing procedure
!
          indxc = 1.D0/(dx(i)+(dx(i+1)+dx(i-1))*0.5D0)
          indzc = 1.D0/(dz(k)+(dz(k+1)+dz(k-1))*0.5D0)

          !ugc = 0.5D0 * (b_e * ug(ijk) + b_w * ug(imjk))
          ugc = 0.5D0 * (ug(ijk) + ug(imjk))

          !wgc = 0.5D0 * (b_t * wg(ijk) + b_b * wg(ijkm))
          wgc = 0.5D0 * (wg(ijk) + wg(ijkm))

          IF (job_type == '3D') THEN
            indyc = 1.D0 / (dy(j)+(dy(j+1)+dy(j-1))*0.5D0)

            !vgc = 0.5D0 * (b_n * vg(ijk) + b_s * vg(ijmk))
            vgc = 0.5D0 * (vg(ijk) + vg(ijmk))
          END IF
!
! ... Pressure term (can be coupled in pressure algorithm)
!
          dpxyz= 0.D0
          dpxyz = dpxyz + dt * indxc * ugc * (p(ijke)-p(ijkw))
          dpxyz = dpxyz + dt * indyc * vgc * (p(ijkn)-p(ijks))
          dpxyz = dpxyz + dt * indzc * wgc * (p(ijkt)-p(ijkb))
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
! ... Compute specific heat for gas (cp, cg) and particles (ck);
! ... update temperature of gas (tg ) and particles (ts)
! ... from caloric Equation of State
!
          CALL caloric_eosg(cp(:,ijk), cg(ijk), tg(ijk), ygc(ijk,:), &
                            sieg(ijk), ijk, info)
          num = num + info

          DO is=1, nsolid
            CALL caloric_eosl(ts(ijk,is),cps(is),ck(is,ijk),sies(ijk,is),ijk) 
          END DO

        END IF

      END DO
!
! ... Check the simultaneous convergence of the caloric eos 
! ... for all procs
!
      CALL parallel_sum_integer( num, 1 )
      IF ( num > 0 ) &
        CALL error('ftem','Error in caloric equation of state', num)
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
        IF(rlk(ijk,is-1)*inrl(is-1) <= 1.D-9) THEN
!        IF(ABS(a(is,is)) <= flim) THEN
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
