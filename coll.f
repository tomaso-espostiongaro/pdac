!----------------------------------------------------------------------
      SUBROUTINE collect
!----------------------------------------------------------------------
! ... This routine collects local variables form computing units to 
! ... store global variables, using PARALLEL_SUM MPI routine
!
      USE eos_gas, ONLY: gc_bulk_density, gc_mass_fraction, gc_molar_fraction
      USE eos_gas, ONLY: rgpgc, ygc, xgc
      USE eos_gas, ONLY: gas_heat_capacity
      USE eos_gas, ONLY: cg
      USE gas_solid_velocity, ONLY: gas_velocity_r, gas_velocity_z, gas_velocity_x, gas_velocity_y
      USE gas_solid_velocity, ONLY: solid_velocity_r, solid_velocity_z, solid_velocity_x, solid_velocity_y
      USE gas_solid_velocity, ONLY: ug, wg, us, ws, vg, vs
      USE gas_solid_density, ONLY: gas_bulk_density, solid_bulk_density, gas_density
      USE gas_solid_density, ONLY: rgp, rlk, rog
      USE gas_solid_temperature, ONLY: gas_enthalpy, gas_temperature
      USE gas_solid_temperature, ONLY: solid_enthalpy, solid_temperature
      USE gas_solid_temperature, ONLY: sieg, sies, tg, ts
      USE grid, ONLY: ncint, myijk
      USE pressure_epsilon, ONLY: gas_pressure , void_fraction
      USE pressure_epsilon, ONLY: p, ep
      USE heat_capacity, ONLY: gc_heat_capacity, solid_heat_capacity
      USE heat_capacity, ONLY: cp, ck
      USE control_flags, ONLY: job_type
      USE indijk_module, ONLY: ip0_jp0_kp0_
!
      USE turbulence, ONLY: scoeff, smag_coeff

      IMPLICIT NONE

      INTEGER :: imesh, ijk
!
      gas_pressure                = 0.D0
      solid_bulk_density          = 0.D0
      gas_enthalpy                = 0.D0
      IF( job_type == '2D' ) THEN
        gas_velocity_r              = 0.D0
      ELSE
        gas_velocity_x              = 0.D0
        gas_velocity_y              = 0.D0
      END IF
      gas_velocity_z              = 0.D0
      solid_enthalpy              = 0.D0
      IF( job_type == '2D' ) THEN
        solid_velocity_r            = 0.D0
      ELSE
        solid_velocity_x            = 0.D0
        solid_velocity_y            = 0.D0
      END IF
      solid_velocity_z            = 0.D0
      gc_mass_fraction            = 0.D0
!
      gas_bulk_density            = 0.D0
      gas_density                 = 0.D0
      void_fraction               = 0.D0
      gas_temperature             = 0.D0
      solid_temperature           = 0.D0
      gc_bulk_density             = 0.D0
      gc_molar_fraction           = 0.D0
!
      solid_heat_capacity         = 0.D0
      gc_heat_capacity            = 0.D0
      gas_heat_capacity           = 0.D0
!
      smag_coeff                  = 0.D0
!
      DO ijk = 1, ncint
        imesh = myijk( ip0_jp0_kp0_, ijk)
!
        gas_pressure(imesh) = p(ijk)
        solid_bulk_density(:,imesh) = rlk(:,ijk)
        gas_enthalpy(imesh) = sieg(ijk)
        IF( job_type == '2D' ) THEN
          gas_velocity_r(imesh) = ug(ijk)
        ELSE
          gas_velocity_x(imesh) = ug(ijk)
          gas_velocity_y(imesh) = vg(ijk)
        END IF
        gas_velocity_z(imesh) = wg(ijk)
        solid_enthalpy(:,imesh) = sies(:,ijk)
        IF( job_type == '2D' ) THEN
          solid_velocity_r(:,imesh) = us(:,ijk)
        ELSE
          solid_velocity_x(:,imesh) = us(:,ijk)
          solid_velocity_y(:,imesh) = vs(:,ijk)
        END IF
        solid_velocity_z(:,imesh) = ws(:,ijk)
        gc_mass_fraction(:, imesh) = ygc(:, ijk)
!
        gas_bulk_density(imesh) = rgp(ijk)
        gas_density(imesh) = rog(ijk)
        void_fraction(imesh) = ep(ijk)
        gas_temperature(imesh) = tg(ijk)
        solid_temperature(:,imesh) = ts(:,ijk)
        gc_bulk_density(:, imesh) = rgpgc(:, ijk)
        gc_molar_fraction(:, imesh) = xgc(:, ijk)
!
        solid_heat_capacity(:,imesh) = ck(:,ijk)
        gc_heat_capacity(:,imesh) = cp(:,ijk)
        gas_heat_capacity(imesh) = cg(ijk)
!
        smag_coeff(imesh) = scoeff(ijk)
        
      END DO
!
      CALL parallel_sum_real(gas_pressure, SIZE(gas_pressure))
      CALL parallel_sum_real(solid_bulk_density, SIZE(solid_bulk_density))
      CALL parallel_sum_real(gas_enthalpy, SIZE(gas_enthalpy) )
      IF( job_type == '2D' ) THEN
        CALL parallel_sum_real(gas_velocity_r, SIZE(gas_velocity_r) )
      ELSE
        CALL parallel_sum_real(gas_velocity_x, SIZE(gas_velocity_r) )
        CALL parallel_sum_real(gas_velocity_y, SIZE(gas_velocity_r) )
      END IF
      CALL parallel_sum_real(gas_velocity_z, SIZE(gas_velocity_z) )
      CALL parallel_sum_real(solid_enthalpy, SIZE(solid_enthalpy) )
      IF( job_type == '2D' ) THEN
        CALL parallel_sum_real(solid_velocity_r, SIZE(solid_velocity_r) )
      ELSE
        CALL parallel_sum_real(solid_velocity_x, SIZE(solid_velocity_r) )
        CALL parallel_sum_real(solid_velocity_y, SIZE(solid_velocity_r) )
      END IF
      CALL parallel_sum_real(solid_velocity_z, SIZE(solid_velocity_z) )
      CALL parallel_sum_real(gc_mass_fraction, SIZE(gc_mass_fraction))
!
      CALL parallel_sum_real(gas_bulk_density, SIZE(gas_bulk_density))
      CALL parallel_sum_real(gas_density, SIZE(gas_density))
      CALL parallel_sum_real(void_fraction, SIZE(void_fraction) )
      CALL parallel_sum_real(gas_temperature, SIZE(gas_temperature) )
      CALL parallel_sum_real(solid_temperature, SIZE(solid_temperature) )
      CALL parallel_sum_real(gc_bulk_density, SIZE(gc_bulk_density))
      CALL parallel_sum_real(gc_molar_fraction, SIZE(gc_molar_fraction))
!
      CALL parallel_sum_real(solid_heat_capacity, SIZE(solid_heat_capacity) )
      CALL parallel_sum_real(gc_heat_capacity, SIZE(gc_heat_capacity) )
      CALL parallel_sum_real(gas_heat_capacity, SIZE(gas_heat_capacity) )
!
      CALL parallel_sum_real(smag_coeff, SIZE(smag_coeff) )
!
      END SUBROUTINE
