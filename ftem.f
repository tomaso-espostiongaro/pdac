!----------------------------------------------------------------------
      MODULE enthalpy_matrix
!----------------------------------------------------------------------
      IMPLICIT NONE

! ... Interphase enthalpy matrix elements
!
      PUBLIC
      REAL*8, DIMENSION(:),   ALLOCATABLE :: bt
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: at

      REAL*8 :: flim
      REAL*8, PRIVATE :: hv
      LOGICAL :: tforce
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
      USE immersed_boundaries, ONLY: faces, immb
      USE specific_heat_module, ONLY: ck, cp
      USE heat_transfer, ONLY: hvs
      USE particles_constants, ONLY: cps
      USE pressure_epsilon, ONLY: ep, p, pn
      USE reactions, ONLY: hrex, irex
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm, ipjk, ijpk, ijkp
      USE set_indexes, ONLY: ijke, ijkn, ijkt, ijkw, ijks, ijkb
      USE tilde_energy, ONLY: rhg, rhs
      USE time_parameters, ONLY: time, dt
      IMPLICIT NONE
!
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
        !
        ! ... Compute the new enthalpy in fluid cells 
        ! ... and immersed boundaries
        IF ( BTEST(flag(ijk),0) ) THEN
          CALL meshinds(ijk,imesh,i,j,k)
          CALL subscr(ijk)
!
          IF (immb == 1) CALL faces(ijk, b_e, b_w, b_t, b_b, b_n, b_s, ivf)
          IF (irex == 2) CALL hrex(ijk,hrexg,hrexs)
!
! ... The pressure terms can be treated implicitly
! ... (coupled with momentum equation into the iterative solver)
! ... This term is not modified (?) by the forcing procedure
!
          indxc = 1.D0/(dx(i)+(dx(i+1)+dx(i-1))*0.5D0)
          indzc = 1.D0/(dz(k)+(dz(k+1)+dz(k-1))*0.5D0)

          ugc = 0.5D0 * (b_e * ug(ijk) + b_w * ug(imjk))
          !ugc = 0.5D0 * (ug(ijk) + ug(imjk))

          wgc = 0.5D0 * (b_t * wg(ijk) + b_b * wg(ijkm))
          !wgc = 0.5D0 * (wg(ijk) + wg(ijkm))

          IF (job_type == '3D') THEN
            indyc = 1.D0 / (dy(j)+(dy(j+1)+dy(j-1))*0.5D0)

            vgc = 0.5D0 * (b_n * vg(ijk) + b_s * vg(ijmk))
            !vgc = 0.5D0 * (vg(ijk) + vg(ijmk))
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
! ... update temperature of gas (tg) and particles (ts)
! ... from caloric Equation of State
!
          CALL caloric_eosg(cp(:,ijk), cg(ijk), tg(ijk), ygc(ijk,:), &
                            sieg(ijk), ijk, info)
          num = num + info

          DO is=1, nsolid
            CALL caloric_eosl(ts(ijk,is),cps(is),ck(is,ijk),sies(ijk,is),ijk) 
          END DO

          ! ... Constrain the temperature variation
          !
          IF (tforce) THEN
                  IF (info == 1 ) THEN
                          CALL temperature_limiter(ijk)
                  END IF
          END IF
          
        END IF
      END DO
!
! ... Check the simultaneous convergence of the caloric eos 
! ... for all procs.
! ... Crash if the number of cells not converged exceeds 1000
!
      CALL parallel_sum_integer( num, 1 )
      IF ( num > 1000 ) &
        CALL error('ftem','Error in caloric equation of state', num)
!
      DEALLOCATE(at)
      DEALLOCATE(bt)
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE temperature_limiter(ijk)

      USE control_flags, ONLY: job_type
      USE dimensions, ONLY: nsolid, ngas
      USE eos_gas, ONLY: ygc, cg
      USE gas_constants, ONLY: gas_type, gmw, rgas, tzero, hzerog, hzeros
      USE gas_solid_temperature, ONLY: tg, ts, sieg, sies
      USE io_files, ONLY: testunit
      USE set_indexes, ONLY: imjk, ijmk, ijkm, ipjk, ijpk, ijkp
      USE specific_heat_module, ONLY: ck, cp, hcapg, hcaps
      USE particles_constants, ONLY: inrl, cps

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ijk
      INTEGER :: ig, is
      REAL*8 :: av_tg
      REAL*8 :: av_ts
      REAL*8 :: hc

      IF (job_type == '3D') THEN
        !
        ! ... Gas Temperature and Enthalpy
          av_tg = tg(imjk) + tg(ijmk) + tg(ijkm) + tg(ipjk) + tg(ijpk) + tg(ijkp)
          av_tg = av_tg / 6.D0
          tg(ijk) = av_tg
          !
          CALL hcapg(cp(:,ijk), tg(ijk))
          hc = 0.D0
          DO ig = 1, ngas
            hc = hc + cp(gas_type(ig),ijk) * ygc(ijk,ig)
          END DO
          cg(ijk) = hc
          sieg(ijk) = (tg(ijk)-tzero) * cg(ijk) + hzerog
          WRITE(testunit,*) 'Limiting gas temperature in cell: ', ijk
        !
        ! ... Solid Temperature and Enthalpy
        DO is = 1, nsolid
            av_ts = ts(imjk,is) + ts(ijmk,is) + ts(ijkm,is) + ts(ipjk,is) + ts(ijpk,is) + ts(ijkp,is)
            av_ts = av_ts / 6.D0
            ts(ijk,is) = av_ts
            !
            CALL hcaps(ck(is,ijk), cps(is), ts(ijk,is))
            sies(ijk,is) = ( ts(ijk,is) - tzero ) * ck(is,ijk) + hzeros
            WRITE(testunit,*) 'Limiting particle ',is,' temperature in cell: ', ijk
        END DO
        !
      ELSE IF (job_type == '2D') THEN
        !
        ! ... Gas Temperature and Enthalpy
          av_tg = tg(imjk) + tg(ijkm) + tg(ipjk) + tg(ijkp)
          av_tg = 0.25D0 * av_tg
          tg(ijk) = av_tg
          !
          CALL hcapg(cp(:,ijk), tg(ijk))
          hc = 0.D0
          DO ig = 1, ngas
            hc = hc + cp(gas_type(ig),ijk) * ygc(ijk,ig)
          END DO
          cg(ijk) = hc
          sieg(ijk) = (tg(ijk)-tzero) * cg(ijk) + hzerog
          WRITE(testunit,*) 'Limiting gas temperature in cell: ', ijk
        !
        ! ... Solid Temperature and Enthalpy
        DO is = 1, nsolid
            av_ts = ts(imjk,is) + ts(ijkm,is) + ts(ipjk,is) + ts(ijkp,is)
            av_ts = 0.25D0 * av_ts 
            ts(ijk,is) = av_ts
            !
            CALL hcaps(ck(is,ijk), cps(is), ts(ijk,is))
            sies(ijk,is) = ( ts(ijk,is) - tzero ) * ck(is,ijk) + hzeros
            WRITE(testunit,*) 'Limiting particle ', is, ' temperature in cell: ', ijk
        END DO
        !
      END IF

      RETURN
      END SUBROUTINE temperature_limiter
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
        !IF(rlk(ijk,is-1)*inrl(is-1) < flim*1.D-3) THEN
        IF(ABS(a(is,is)) < flim) THEN
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
