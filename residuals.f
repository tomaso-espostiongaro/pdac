!----------------------------------------------------------------------
      MODULE check_residuals
!
! ... Check the residuals of mass in the whole domain
!
! ... Authors: Tomaso Esposti Ongaro
! ... Date   : 
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
      USE dimensions, ONLY: ngas, nsolid
      USE domain_mapping, ONLY: meshinds
      USE eos_gas, ONLY: ygc
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE grid, ONLY: dx, dy, dz, r, rb, flag
      USE grid, ONLY: inlet_cell, vent_cell
      USE immersed_boundaries, ONLY: immb, faces
      USE io_files, ONLY: checkunit
      USE time_parameters, ONLY: dt
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nswp
      INTEGER :: ijk, i, j, k, imesh
      INTEGER :: ig, is
      INTEGER :: b_e, b_w, b_t, b_b, b_n, b_s
      REAL*8 :: ivf
      REAL*8 :: volume, sx, sy, sz, flux
      REAL*8 :: res_g, mfr
      REAL*8, ALLOCATABLE :: res_s(:), res_gc(:)
      !
      ! ... Initialize the cell fractions for immersed boundaries
      !
      b_e = 1; b_w = 1; b_t = 1; b_b = 1; b_n = 1; b_s = 1; ivf = 1.D0

      ALLOCATE(res_s(nsolid))
      ALLOCATE(res_gc(ngas))

      res_g  = 0.D0
      res_s  = 0.D0
      res_gc = 0.D0

      DO ijk = 1, ncint
        CALL meshinds(ijk,imesh,i,j,k)

        IF ( BTEST(flag(ijk),0) ) THEN
          !
          ! ... Compute the volumes partially filled by the
          ! ... topography
          !
          IF (immb == 1) CALL faces(ijk, b_e, b_w, b_t, b_b, b_n, b_s, ivf)
  
          IF (job_type == '2D') THEN
            volume = r(i) * dx(i) * dz(k) / ivf
          ELSE IF (job_type == '3D') THEN
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
          ! ... Compute the mass entered since the beginning
          !
          IF (job_type == '2D') THEN
            sx = dz(k)*rb(i)
            sy = 0.D0
            sz = dx(i)*r(i)
          ELSE IF (job_type == '3D') THEN
            sx = dz(k)*dy(j)
            sy = dx(i)*dz(k)
            sz = dy(j)*dx(i)
          END IF
          
          flux = ug(ijk)*sx + wg(ijk) * sz
          IF (job_type == '3D') flux = flux + vg(ijk) * sy
          flux = rgp(ijk) * flux
          res_g = res_g - flux * nswp * dt

          DO is = 1, nsolid
            flux = us(ijk,is)*sx + ws(ijk,is) * sz
            IF (job_type == '3D') flux = flux + vs(ijk,is) * sy
            flux = rlk(ijk,is) * flux
            res_s(is) = res_s(is) - flux * nswp * dt
          END DO
          
          DO ig = 1, ngas
            flux = ug(ijk)*sx + wg(ijk) * sz
            IF (job_type == '3D') flux = flux + vg(ijk) * sy
            flux = rgp(ijk) * ygc(ijk,ig) * flux
            res_gc(ig) = res_gc(ig) - flux * nswp * dt
          END DO

        END IF

      END DO

      CALL parallel_sum_real(res_g, 1)
      CALL parallel_sum_real(res_s, nsolid)
      CALL parallel_sum_real(res_gc,ngas)

      CALL compute_mass_flow_rate(mfr)

      IF (mpime == root) THEN
        WRITE(checkunit,55) nswp, res_g, (res_s(is),  is=1,nsolid), &
                                  (res_gc(ig), ig=1, ngas), mfr
      END IF

 55   FORMAT(I8,15(G30.20E3))

      DEALLOCATE(res_s)
      DEALLOCATE(res_gc)

      RETURN
      END SUBROUTINE print_mass_residuals
!----------------------------------------------------------------------
      SUBROUTINE compute_mass_flow_rate(mfr)

      USE control_flags, ONLY: job_type
      USE dimensions, ONLY: ngas, nsolid
      USE domain_mapping, ONLY: meshinds
      USE eos_gas, ONLY: ygc
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: tg, ts
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE grid, ONLY: dx, dy, dz, r, rb, flag
      USE grid, ONLY: inlet_cell, vent_cell
      USE pressure_epsilon, ONLY: p, ep
      USE time_parameters, ONLY: dt
      IMPLICIT NONE

      INTEGER :: ijk, i, j, k, imesh
      INTEGER :: ig, is
      REAL*8 :: volume, sx, sy, sz, flux
      REAL*8, INTENT(OUT) :: mfr

      mfr  = 0.D0

      DO ijk = 1, ncint
        CALL meshinds(ijk,imesh,i,j,k)

        IF (flag(ijk) == vent_cell) THEN
          !
          ! ... Compute the mass entered since the beginning
          !
          IF (job_type == '2D') THEN
            sx = dz(k)*rb(i)
            sy = 0.D0
            sz = dx(i)*r(i)
          ELSE IF (job_type == '3D') THEN
            sx = dz(k)*dy(j)
            sy = dx(i)*dz(k)
            sz = dy(j)*dx(i)
          END IF
          
          flux = rgp(ijk) * wg(ijk) * sz
          mfr = mfr + flux

          DO is = 1, nsolid
            flux = rlk(ijk,is) * ws(ijk,is) * sz
            mfr = mfr + flux
          END DO
          
        END IF

      END DO

      CALL parallel_sum_real(mfr, 1)

      RETURN
      END SUBROUTINE compute_mass_flow_rate
!----------------------------------------------------------------------
      END MODULE check_residuals
!----------------------------------------------------------------------
