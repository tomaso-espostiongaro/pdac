!-----------------------------------------------------------------------
      SUBROUTINE distribute
!-----------------------------------------------------------------------
! ... This routine distribute global initial variables among computing 
! ... units to store local variables (inverse of `collect' routine)
!
      USE eos_gas, ONLY: gc_bulk_density, gc_mass_fraction, gc_molar_fraction
      USE eos_gas, ONLY: rgpgc, ygc, xgc
      USE eos_gas, ONLY: gas_heat_capacity
      USE eos_gas, ONLY: cg
      USE gas_solid_velocity, ONLY: gas_velocity_r, gas_velocity_z, gas_velocity_x, gas_velocity_y
      USE gas_solid_velocity, ONLY: solid_velocity_r, solid_velocity_z, solid_velocity_x, solid_velocity_y
      USE gas_solid_velocity, ONLY: ug, vg, wg, uk, vk, wk
      USE gas_solid_density, ONLY: gas_bulk_density, solid_bulk_density, gas_density
      USE gas_solid_density, ONLY: rgp, rlk, rog
      USE gas_solid_temperature, ONLY: gas_enthalpy, gas_temperature
      USE gas_solid_temperature, ONLY: solid_enthalpy, solid_temperature
      USE gas_solid_temperature, ONLY: tg, tk, sieg, siek
      USE gas_solid_viscosity, ONLY: mus, particle_viscosity
      USE grid, ONLY: nij_l, myij, data_exchange
      USE pressure_epsilon, ONLY: gas_pressure, void_fraction
      USE pressure_epsilon, ONLY: p, ep
      USE heat_capacity, ONLY: gc_heat_capacity, solid_heat_capacity
      USE heat_capacity, ONLY: cp, ck
      USE turbulence, ONLY: smag, smag_factor
      USE control_flags, ONLY: job_type
!
      INTEGER :: ij, ij_l
      IMPLICIT NONE
!       
! ... distribute arrays among processors
!
        p = 0.D0
        rlk = 0.D0
        sieg = 0.D0
        ug = 0.D0
        vg = 0.D0
        siek = 0.D0
        uk = 0.D0
        vk = 0.D0
        ygc = 0.D0

        IF( job_type == '3D' ) THEN
          wg = 0.D0
          wk = 0.D0
        END IF
!
        rgp = 0.D0
        rog = 0.D0
        ep = 0.D0
        tg = 0.D0
        tk = 0.D0
        rgpgc = 0.D0
        xgc = 0.D0
!
! ... main constitutive parameters
!
        cp = 0.D0
        ck = 0.D0
        cg = 0.D0
        smag = 0.D0
        mus = 0.D0
!        
      DO ij_l = 1, nij_l
        ij = myij(0, 0, ij_l)
!
        p(ij_l) = gas_pressure(ij)
        rlk(:,ij_l) = solid_bulk_density(:,ij)
        sieg(ij_l) = gas_enthalpy(ij)

        IF( job_type == '2D' ) THEN
          ug(ij_l) = gas_velocity_r(ij)
          vg(ij_l) = gas_velocity_z(ij)
        ELSE
          ug(ij_l) = gas_velocity_x(ij)
          vg(ij_l) = gas_velocity_y(ij)
          wg(ij_l) = gas_velocity_z(ij)
        END IF

        siek(:,ij_l) = solid_enthalpy(:,ij)

        IF( job_type == '2D' ) THEN
          uk(:,ij_l) = solid_velocity_r(:,ij)
          vk(:,ij_l) = solid_velocity_z(:,ij)
        ELSE
          uk(:,ij_l) = solid_velocity_x(:,ij)
          vk(:,ij_l) = solid_velocity_y(:,ij)
          wk(:,ij_l) = solid_velocity_z(:,ij)
        END IF

        ygc(:, ij_l) = gc_mass_fraction(:, ij)
!
        rgp(ij_l) = gas_bulk_density(ij)
        rog(ij_l) = gas_density(ij)
        ep(ij_l) = void_fraction(ij)
        tg(ij_l) = gas_temperature(ij)
        tk(:,ij_l) = solid_temperature(:,ij)
        rgpgc(:, ij_l) = gc_bulk_density(:, ij)
        xgc(:, ij_l) = gc_molar_fraction(:, ij)
!
        ck(:,ij_l) = solid_heat_capacity(:,ij)
        cp(:,ij_l) =gc_heat_capacity(:,ij)
        cg(ij_l) = gas_heat_capacity(ij)
        smag(ij_l) = smag_factor(ij)
        mus(:,ij_l) = particle_viscosity(:,ij)
      END DO
!
      END SUBROUTINE
