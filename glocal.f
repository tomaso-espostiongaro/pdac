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
      USE specific_heat, ONLY: gc_specific_heat, solid_specific_heat
      USE specific_heat, ONLY: cp, ck

      USE turbulence_model, ONLY: smag, smag_factor
      USE turbulence_model, ONLY: scoeff, smag_coeff
!
      IMPLICIT NONE
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE collect
!----------------------------------------------------------------------
!
! ... This routine collects local variables form processors to 
! ... store global variables, using PARALLEL_SUM MPI routine
!
      USE control_flags, ONLY: job_type
      USE grid, ONLY: ncint, myijk
      USE indijk_module, ONLY: ip0_jp0_kp0_
      IMPLICIT NONE

      INTEGER :: ijk, imesh
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
      solid_specific_heat         = 0.D0
      gc_specific_heat            = 0.D0
      gas_specific_heat           = 0.D0
!
      smag_coeff                  = 0.D0
!
      DO ijk = 1, ncint
        imesh = myijk( ip0_jp0_kp0_, ijk)
!
        gas_pressure(imesh) = p(ijk)
        solid_bulk_density(imesh,:) = rlk(ijk,:)
        gas_enthalpy(imesh) = sieg(ijk)
        IF( job_type == '2D' ) THEN
          gas_velocity_r(imesh) = ug(ijk)
        ELSE
          gas_velocity_x(imesh) = ug(ijk)
          gas_velocity_y(imesh) = vg(ijk)
        END IF
        gas_velocity_z(imesh) = wg(ijk)
        solid_enthalpy(imesh,:) = sies(ijk,:)
        IF( job_type == '2D' ) THEN
          solid_velocity_r(imesh,:) = us(ijk,:)
        ELSE
          solid_velocity_x(imesh,:) = us(ijk,:)
          solid_velocity_y(imesh,:) = vs(ijk,:)
        END IF
        solid_velocity_z(imesh,:) = ws(ijk,:)
        gc_mass_fraction(:, imesh) = ygc(:, ijk)
!
        gas_bulk_density(imesh) = rgp(ijk)
        gas_density(imesh) = rog(ijk)
        void_fraction(imesh) = ep(ijk)
        gas_temperature(imesh) = tg(ijk)
        solid_temperature(imesh,:) = ts(ijk,:)
        gc_bulk_density(imesh,:) = rgpgc(ijk,:)
        gc_molar_fraction(:, imesh) = xgc(:,ijk)
!
        solid_specific_heat(:,imesh) = ck(:,ijk)
        gc_specific_heat(:,imesh) = cp(:,ijk)
        gas_specific_heat(imesh) = cg(ijk)
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
        CALL parallel_sum_real(gas_velocity_x, SIZE(gas_velocity_x) )
        CALL parallel_sum_real(gas_velocity_y, SIZE(gas_velocity_y) )
      END IF
      CALL parallel_sum_real(gas_velocity_z, SIZE(gas_velocity_z) )
      CALL parallel_sum_real(solid_enthalpy, SIZE(solid_enthalpy) )
      IF( job_type == '2D' ) THEN
        CALL parallel_sum_real(solid_velocity_r, SIZE(solid_velocity_r) )
      ELSE
        CALL parallel_sum_real(solid_velocity_x, SIZE(solid_velocity_x) )
        CALL parallel_sum_real(solid_velocity_y, SIZE(solid_velocity_y) )
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
      CALL parallel_sum_real(solid_specific_heat, SIZE(solid_specific_heat) )
      CALL parallel_sum_real(gc_specific_heat, SIZE(gc_specific_heat) )
      CALL parallel_sum_real(gas_specific_heat, SIZE(gas_specific_heat) )
!
      CALL parallel_sum_real(smag_coeff, SIZE(smag_coeff) )
!
      END SUBROUTINE collect
!-----------------------------------------------------------------------
      SUBROUTINE distribute

! ... This routine distribute global initial variables among processors
! ... to store local variables (inverse of `collect' routine)
!
      USE control_flags, ONLY: job_type
      USE grid, ONLY: ncint, myijk
      USE indijk_module, ONLY: ip0_jp0_kp0_
!
      IMPLICIT NONE
!
      INTEGER :: ijk, imesh
!       
! ... distribute arrays among processors
!
        p = 0.D0
        rlk = 0.D0
        sieg = 0.D0
        ug = 0.D0
        wg = 0.D0
        sies = 0.D0
        us = 0.D0
        ws = 0.D0
        ygc = 0.D0

        IF( job_type == '3D' ) THEN
          vg = 0.D0
          vs = 0.D0
        END IF
!
        rgp = 0.D0
        rog = 0.D0
        ep = 0.D0
        tg = 0.D0
        ts = 0.D0
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
      DO ijk = 1, ncint
        imesh = myijk( ip0_jp0_kp0_, ijk)
!
        p(ijk) = gas_pressure(imesh)
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
!
        rgp(ijk) = gas_bulk_density(imesh)
        rog(ijk) = gas_density(imesh)
        ep(ijk) = void_fraction(imesh)
        tg(ijk) = gas_temperature(imesh)
        ts(ijk,:) = solid_temperature(imesh,:)
        rgpgc(ijk,:) = gc_bulk_density(imesh,:)
        xgc(:,ijk) = gc_molar_fraction(:, imesh)
!
        ck(:,ijk) = solid_specific_heat(:,imesh)
        cp(:,ijk) = gc_specific_heat(:,imesh)
        cg(ijk) = gas_specific_heat(imesh)
        smag(ijk) = smag_factor(imesh)
        mus(ijk,:) = particle_viscosity(imesh,:)
      END DO
!
      END SUBROUTINE distribute
!-----------------------------------------------------------------------
      END MODULE glocal_arrays
!-----------------------------------------------------------------------
