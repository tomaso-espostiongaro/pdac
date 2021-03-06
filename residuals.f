!----------------------------------------------------------------------
      MODULE check_residuals
!
! ... Integrate the gas and solid densities over the whole domain,
! ... including a correction for immersed-boundary cells. 
! ... Compute the istantaneous mass flow-rate from a vent.
!----------------------------------------------------------------------
      USE domain_mapping, ONLY: ncint
      USE parallel, ONLY: root, mpime
      IMPLICIT NONE
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE print_mass_residuals(nswp)

      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: ngas, nsolid
      USE domain_mapping, ONLY: meshinds
      USE eos_gas, ONLY: ygc
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE grid, ONLY: dx, dy, dz, r, rb, flag
      USE grid, ONLY: inlet_cell, vent_cell
      USE immersed_boundaries, ONLY: faces
      USE io_files, ONLY: checkunit
      USE time_parameters, ONLY: dt
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nswp
      INTEGER :: ijk, i, j, k, imesh
      INTEGER :: ig, is
      INTEGER :: b_e, b_w, b_t, b_b, b_n, b_s
      REAL*8 :: ivf
      REAL*8 :: volume, sx, sy, sz, flux
      REAL*8 :: res_g, mfr, mgd, mxv, mrd
      REAL*8, ALLOCATABLE :: res_s(:), res_gc(:)
      REAL*8, ALLOCATABLE :: msd(:)
!
      ALLOCATE(msd(nsolid))
      ALLOCATE(res_s(nsolid))
      ALLOCATE(res_gc(ngas))

      res_g  = 0.D0
      res_s  = 0.D0
      res_gc = 0.D0

      DO ijk = 1, ncint
        CALL meshinds(ijk,imesh,i,j,k)

        ! ... On computational fluid cells, compute the density and the volume
        ! ... to obtain the total mass
        !
        IF ( BTEST(flag(ijk),0) ) THEN
          !
          ! ... Compute the volumes partially filled by the
          ! ... topography (IVF)
          !
          CALL faces(ijk, b_e, b_w, b_t, b_b, b_n, b_s, ivf)
  
          IF (job_type == JOB_TYPE_2D) THEN
            volume = r(i) * dx(i) * dz(k) / ivf
          ELSE IF (job_type == JOB_TYPE_3D) THEN
            volume = dx(i) * dy(j) * dz(k) / ivf
          END IF
  
          ! ... Compute the total mass in the physical domain
          !
          res_g     = res_g + rgp(ijk) * volume

          DO is = 1, nsolid
            res_s(is)  = res_s(is) + rlk(ijk,is) * volume
          END DO

          DO ig = 1, ngas
            res_gc(ig) = res_gc(ig) + rgp(ijk) * ygc(ijk,ig) * volume
          END DO

        ELSE IF (flag(ijk) == inlet_cell .OR. flag(ijk) == vent_cell) THEN
          !
          ! ... Compute the mass entered since the beginning and subtract from 
          ! ... the total mass
          !
          IF (job_type == JOB_TYPE_2D) THEN
            sx = dz(k)*rb(i)
            sy = 0.D0
            sz = dx(i)*r(i)
          ELSE IF (job_type == JOB_TYPE_3D) THEN
            sx = dz(k)*dy(j)
            sy = dx(i)*dz(k)
            sz = dy(j)*dx(i)
          END IF
          
          flux = ug(ijk)*sx + wg(ijk) * sz
          IF (job_type == JOB_TYPE_3D) flux = flux + vg(ijk) * sy
          flux = rgp(ijk) * flux
          res_g = res_g - flux * nswp * dt

          DO is = 1, nsolid
            flux = us(ijk,is)*sx + ws(ijk,is) * sz
            IF (job_type == JOB_TYPE_3D) flux = flux + vs(ijk,is) * sy
            flux = rlk(ijk,is) * flux
            res_s(is) = res_s(is) - flux * nswp * dt
          END DO
          
          DO ig = 1, ngas
            flux = ug(ijk)*sx + wg(ijk) * sz
            IF (job_type == JOB_TYPE_3D) flux = flux + vg(ijk) * sy
            flux = rgp(ijk) * ygc(ijk,ig) * flux
            res_gc(ig) = res_gc(ig) - flux * nswp * dt
          END DO

        END IF

      END DO

      CALL parallel_sum_real(res_g, 1)
      CALL parallel_sum_real(res_s, nsolid)
      CALL parallel_sum_real(res_gc,ngas)
!
      CALL compute_mass_flow_rate(mfr, mgd, msd, mxv, mrd)
!
      IF (mpime == root) THEN
        WRITE(checkunit,55) nswp, res_g, (res_s(is),  is=1,nsolid), &
                                  (res_gc(ig), ig=1, ngas), mfr
      END IF

 55   FORMAT(I8,15(G30.20E3))

      DEALLOCATE(msd)
      DEALLOCATE(res_s)
      DEALLOCATE(res_gc)

      RETURN
      END SUBROUTINE print_mass_residuals
!----------------------------------------------------------------------
      SUBROUTINE print_mass_flow_rate(mfr)
      USE control_flags, ONLY: lpr
      USE dimensions
      USE io_files, ONLY: logunit, ventunit, ventfile
      USE vent_conditions, ONLY: ivent
      IMPLICIT NONE
      REAL*8 :: mgd, mxv, mrd
      REAL*8, INTENT(OUT) :: mfr
      REAL*8, ALLOCATABLE :: msd(:)
!
      ALLOCATE(msd(nsolid))
!
      CALL compute_mass_flow_rate(mfr, mgd, msd, mxv, mrd)
!
      IF (mpime == root) THEN
        IF (mfr > 0.D0) THEN
          WRITE(ventunit,300) 'Mass flow rate          : ', mfr
          WRITE(ventunit,400) 'Gas Density at vent     : ', mgd
          WRITE(ventunit,400) 'Solid Density at vent   : ', msd
          WRITE(ventunit,400) 'Mixture Velocity at vent: ', mxv 
          WRITE(ventunit,400) 'Averaged vent radius    : ', mrd 
          WRITE(ventunit,400) ''
        END IF
      END IF
!
 300    FORMAT(A26,10(F15.5))
 400    FORMAT(A26,10(F12.5))
!
      DEALLOCATE(msd)
!
      RETURN
      END SUBROUTINE print_mass_flow_rate
!----------------------------------------------------------------------
      SUBROUTINE compute_mass_flow_rate(mfr, mgd, msd, mxv, mrd)

      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: ngas, nsolid
      USE domain_mapping, ONLY: meshinds
      USE eos_gas, ONLY: ygc, xgc
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: tg, ts
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE grid, ONLY: dx, dy, dz, r, rb, flag, itc
      USE grid, ONLY: inlet_cell, vent_cell
      USE io_files, ONLY: testunit
      USE pressure_epsilon, ONLY: p, ep
      USE time_parameters, ONLY: dt
      IMPLICIT NONE

      INTEGER :: ijk, i, j, k, imesh
      INTEGER :: ig, is
      REAL*8 :: volume, sx, sy, sz, surface
      REAL*8 :: flux, pi, twopi, mixd
      REAL*8, INTENT(OUT) :: mfr, mgd, mxv, mrd
      REAL*8, INTENT(OUT), DIMENSION(:) :: msd
!
      pi = 4.D0 * ATAN(1.D0)
      twopi = 2.D0 * pi
!            
      surface = 0.D0
      mfr  = 0.D0
      mgd = 0.D0
      msd = 0.D0
      mxv = 0.D0
!
      IF (job_type == JOB_TYPE_2D .AND. itc == 0) RETURN
!
      DO ijk = 1, ncint
        CALL meshinds(ijk,imesh,i,j,k)

        IF (flag(ijk) == vent_cell .OR. flag(ijk) == inlet_cell) THEN
          !
          ! ... Compute the mass entered since the beginning
          !
          IF (job_type == JOB_TYPE_2D .AND. itc==1) THEN
            sx = dz(k)*rb(i)
            sy = 0.D0
            sz = pi*dx(i)*(2.D0*rb(i) - dx(i))
          ELSE IF (job_type == JOB_TYPE_3D) THEN
            sx = dz(k)*dy(j)
            sy = dx(i)*dz(k)
            sz = dy(j)*dx(i)
          END IF
          
          flux = rgp(ijk) * wg(ijk) * sz
          mgd = mgd + rgp(ijk) * sz

          DO is = 1, nsolid
            flux = flux + rlk(ijk,is) * ws(ijk,is) * sz
            msd(is) = msd(is) + rlk(ijk,is) * sz
          END DO
          mfr = mfr + flux
          IF (flux > 0.D0) surface = surface + sz

        END IF

      END DO

      CALL parallel_sum_real(surface, 1)
      CALL parallel_sum_real(mfr, 1)
      CALL parallel_sum_real(mgd, 1)
      CALL parallel_sum_real(msd, nsolid)
      
      IF (surface /= 0.D0) THEN
        mgd = mgd / surface
        msd = msd / surface
        mixd = mgd + SUM(msd)
        mxv = mfr / mixd / surface
        mrd = DSQRT(surface / pi)
      END IF

      RETURN
      END SUBROUTINE compute_mass_flow_rate
!----------------------------------------------------------------------
      END MODULE check_residuals
!----------------------------------------------------------------------
