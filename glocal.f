!----------------------------------------------------------------------
      MODULE glocal_arrays
!
      USE eos_gas, ONLY: gc_bulk_density, gc_mass_fraction, gc_molar_fraction
      USE eos_gas, ONLY: rgpgc, ygc, xgc
      USE eos_gas, ONLY: gas_specific_heat
      USE eos_gas, ONLY: cg
      USE gas_solid_velocity, ONLY: gas_velocity_r, gas_velocity_z,      &
                                    gas_velocity_x, gas_velocity_y
      USE gas_solid_velocity, ONLY: solid_velocity_r, solid_velocity_z,  &
                                    solid_velocity_x, solid_velocity_y
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_viscosity, ONLY: mus, particle_viscosity
      USE gas_solid_density, ONLY: gas_bulk_density, solid_bulk_density, &
                                   gas_density
      USE gas_solid_density, ONLY: rgp, rlk, rog
      USE gas_solid_temperature, ONLY: gas_enthalpy, gas_temperature
      USE gas_solid_temperature, ONLY: solid_enthalpy, solid_temperature
      USE gas_solid_temperature, ONLY: sieg, sies, tg, ts
      USE pressure_epsilon, ONLY: gas_pressure , void_fraction
      USE pressure_epsilon, ONLY: p, ep
      USE specific_heat_module, ONLY: gc_specific_heat, solid_specific_heat
      USE specific_heat_module, ONLY: cp, ck

      USE turbulence_model, ONLY: smag, smag_factor
      USE turbulence_model, ONLY: scoeff, smag_coeff
!
      IMPLICIT NONE
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------

      SUBROUTINE distribute( itd )

! ... This routine distribute global initial variables among processors
! ... to store local variables (inverse of `collect' routine)
!
      USE control_flags, ONLY: job_type
      USE grid, ONLY: ncint, myijk
      USE indijk_module, ONLY: ip0_jp0_kp0_
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: itd
      INTEGER :: ijk, imesh
!       
! ... distribute arrays among processors
!

!
!
! ... main constitutive parameters
!
      smag = 0.D0
      mus = 0.D0

      IF( itd <= 1 ) THEN
        p = 0.0d0
        tg = 0.D0
        rlk = 0.D0
        sieg = 0.D0
        ug = 0.D0
        wg = 0.D0
        IF( job_type == '3D' ) THEN
          vg = 0.D0
        END IF
        sies = 0.D0
        us = 0.D0
        ws = 0.D0
        IF( job_type == '3D' ) THEN
          vs = 0.D0
        END IF
        ygc = 0.D0
        rgp = 0.D0
        rog = 0.D0
        ep = 0.D0
        ts = 0.D0
        rgpgc = 0.D0
        xgc = 0.D0
        cg = 0.D0
        cp = 0.D0
        ck = 0.D0
        DO ijk = 1, ncint
          imesh = myijk( ip0_jp0_kp0_, ijk)
          p(ijk) = gas_pressure( imesh )
          tg(ijk) = gas_temperature( imesh )
          rlk(ijk,:) = solid_bulk_density(imesh,:)
          sieg(ijk) = gas_enthalpy(imesh)
          IF( job_type == '2D' ) THEN
            ug(ijk) = gas_velocity_r(imesh)
          ELSE
            ug(ijk) = gas_velocity_x(imesh)
            vg(ijk) = gas_velocity_y(imesh)
          END IF
          wg(ijk) = gas_velocity_z(imesh)
          sies(ijk,:) = solid_enthalpy(imesh,:)
          IF( job_type == '2D' ) THEN
            us(ijk,:) = solid_velocity_r(imesh,:)
          ELSE
            us(ijk,:) = solid_velocity_x(imesh,:)
            vs(ijk,:) = solid_velocity_y(imesh,:)
          END IF
          ws(ijk,:) = solid_velocity_z(imesh,:)
          ygc(:, ijk) = gc_mass_fraction(:, imesh)
          rgp(ijk) = gas_bulk_density(imesh)
          rog(ijk) = gas_density(imesh)
          ep(ijk) = void_fraction(imesh)
          ts(ijk,:) = solid_temperature(imesh,:)
          rgpgc(ijk,:) = gc_bulk_density(imesh,:)
          xgc(:,ijk) = gc_molar_fraction(:, imesh)
          cg(ijk) = gas_specific_heat(imesh)
          ck(:,ijk) = solid_specific_heat(:,imesh)
          cp(:,ijk) = gc_specific_heat(:,imesh)
        END DO
      END IF
!
!        
      DO ijk = 1, ncint
        imesh = myijk( ip0_jp0_kp0_, ijk)
        smag(ijk)  = smag_factor(imesh)
        mus(ijk,:) = particle_viscosity(imesh,:)
      END DO
!
      END SUBROUTINE distribute
!-----------------------------------------------------------------------
      END MODULE glocal_arrays
!-----------------------------------------------------------------------
