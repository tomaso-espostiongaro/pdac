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
      USE gas_solid_velocity, ONLY: gas_velocity_r, gas_velocity_z
      USE gas_solid_velocity, ONLY: solid_velocity_r, solid_velocity_z
      USE gas_solid_velocity, ONLY: ug, vg, uk, vk
      USE gas_solid_density, ONLY: gas_bulk_density, solid_bulk_density, gas_density
      USE gas_solid_density, ONLY: rgp, rlk, rog
      USE gas_solid_temperature, ONLY: gas_enthalpy, gas_temperature
      USE gas_solid_temperature, ONLY: solid_enthalpy, solid_temperature
      USE gas_solid_temperature, ONLY: sieg, siek, tg, tk
      USE grid, ONLY: nij_l, myij
      USE pressure_epsilon, ONLY: gas_pressure, void_fraction
      USE pressure_epsilon, ONLY: p, ep
      USE heat_capacity, ONLY: gc_heat_capacity, solid_heat_capacity
      USE heat_capacity, ONLY: cp, ck
!
      gas_pressure                = 0.D0
      solid_bulk_density          = 0.D0
      gas_enthalpy                = 0.D0
      gas_velocity_r              = 0.D0
      gas_velocity_z              = 0.D0
      solid_enthalpy              = 0.D0
      solid_velocity_r            = 0.D0
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
      DO ij_l = 1, nij_l
        ij = myij(0,0,ij_l)
        gas_pressure(ij) = p(ij_l)
        solid_bulk_density(:,ij) = rlk(:,ij_l)
        gas_enthalpy(ij) = sieg(ij_l)
        gas_velocity_r(ij) = ug(ij_l)
        gas_velocity_z(ij) = vg(ij_l)
        solid_enthalpy(:,ij) = siek(:,ij_l)
        solid_velocity_r(:,ij) = uk(:,ij_l)
        solid_velocity_z(:,ij) = vk(:,ij_l)
        gc_mass_fraction(:, ij) = ygc(:, ij_l)
!
        gas_bulk_density(ij) = rgp(ij_l)
        gas_density(ij) = rog(ij_l)
        void_fraction(ij) = ep(ij_l)
        gas_temperature(ij) = tg(ij_l)
        solid_temperature(:,ij) = tk(:,ij_l)
        gc_bulk_density(:, ij) = rgpgc(:, ij_l)
        gc_molar_fraction(:, ij) = xgc(:, ij_l)
!
        solid_heat_capacity(:,ij) = ck(:,ij_l)
        gc_heat_capacity(:,ij) = cp(:,ij_l)
        gas_heat_capacity(ij) = cg(ij_l)
        
      END DO
!
      CALL parallel_sum_real(gas_pressure, SIZE(gas_pressure))
      CALL parallel_sum_real(solid_bulk_density, SIZE(solid_bulk_density))
      CALL parallel_sum_real(gas_enthalpy, SIZE(gas_enthalpy) )
      CALL parallel_sum_real(gas_velocity_r, SIZE(gas_velocity_r) )
      CALL parallel_sum_real(gas_velocity_z, SIZE(gas_velocity_z) )
      CALL parallel_sum_real(solid_enthalpy, SIZE(solid_enthalpy) )
      CALL parallel_sum_real(solid_velocity_r, SIZE(solid_velocity_r) )
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
      END SUBROUTINE
