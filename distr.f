!-----------------------------------------------------------------------
      SUBROUTINE distribute
!-----------------------------------------------------------------------
! ... This routine distribute global initial variables among computing 
! ... units to store local variables (inverse of `collect' routine)
!
      USE eos_gas, ONLY: rgpgc_g, rgpgcn_g, ygc_g, xgc_g
      USE eos_gas, ONLY: rgpgc, rgpgcn, ygc, xgc
      USE eos_gas, ONLY: cg_g
      USE eos_gas, ONLY: cg
      USE gas_solid_velocity, ONLY: ug_g, vg_g, uk_g, vk_g
      USE gas_solid_velocity, ONLY: ug, vg, uk, vk
      USE gas_solid_density, ONLY: rgp_g, rgpn_g, rlk_g, rlkn_g, rog_g
      USE gas_solid_density, ONLY: rgp, rgpn, rlk, rlkn, rog
      USE gas_solid_temperature, ONLY: siegn_g, siekn_g, sieg_g, siek_g, tg_g, tk_g
      USE gas_solid_temperature, ONLY: siegn, siekn, sieg, siek, tg, tk
      USE grid, ONLY: nij_l, myij, data_exchange
      USE pressure_epsilon, ONLY: p_g, pn_g, ep_g
      USE pressure_epsilon, ONLY: p, pn, ep
      USE th_capacity, ONLY: cp_g, ck_g
      USE th_capacity, ONLY: cp, ck
      USE turbulence, ONLY: smagl_g, smagl
      USE gas_solid_viscosity, ONLY: mug_g, kapg_g
      USE gas_solid_viscosity, ONLY: mug, kapg
      USE turbulence, ONLY: mus, mus_g
!       
! ... distribute the arrays of the main physical quantities among processors
!
        rgp = 0.0
        rgpn = 0.0
        rog = 0.0
        rlk = 0.0
        rlkn = 0.0
        ep = 0.0
        p = 0.0
        pn = 0.0
        sieg = 0.0
        siegn = 0.0
        tg = 0.0
        siek = 0.0
        siekn = 0.0
        tk = 0.0
        rgpgcn = 0.0
        rgpgc = 0.0
        ygc = 0.0
        xgc = 0.0
        ug = 0.0
        vg = 0.0
        uk = 0.0
        vk = 0.0
!
! ... main constitutive paramenters
!
        cp = 0.0
        ck = 0.0
        cg = 0.0
        mug = 0.0
        kapg = 0.0
        mus = 0.0
        smagl = 0.0
!        
      DO ij_l = 1, nij_l
        ij = myij(0, 0, ij_l)
!
        rgp(ij_l) = rgp_g(ij)
        rgpn(ij_l) = rgpn_g(ij)
        rog(ij_l) = rog_g(ij)
        rlk(:,ij_l) = rlk_g(:,ij)
        rlkn(:,ij_l) = rlkn_g(:,ij)
        ep(ij_l) = ep_g(ij)
        p(ij_l) = p_g(ij)
        pn(ij_l) = pn_g(ij)
        sieg(ij_l) = sieg_g(ij)
        siegn(ij_l) = siegn_g(ij)
        tg(ij_l) = tg_g(ij)
        siek(:,ij_l) = siek_g(:,ij)
        siekn(:,ij_l) = siekn_g(:,ij)
        tk(:,ij_l) = tk_g(:,ij)
        rgpgcn(:, ij_l) = rgpgcn_g(:, ij)
        rgpgc(:, ij_l) = rgpgc_g(:, ij)
        ygc(:, ij_l) = ygc_g(:, ij)
        xgc(:, ij_l) = xgc_g(:, ij)
        ug(ij_l) = ug_g(ij)
        vg(ij_l) = vg_g(ij)
        uk(:,ij_l) = uk_g(:,ij)
        vk(:,ij_l) = vk_g(:,ij)
!
        cp(:,ij_l) = cp_g(:,ij)
        ck(:,ij_l) = ck_g(:,ij)
        cg(ij_l) = cg_g(ij)
        mug(ij_l) = mug_g(ij)
        kapg(ij_l) = kapg_g(ij)
        mus(:,ij_l) = mus_g(:,ij)
        smagl(ij_l) = smagl_g(ij)
      END DO
!
      END SUBROUTINE
