!----------------------------------------------------------------------
      SUBROUTINE collect
!----------------------------------------------------------------------
! ... This routine collects local variables form computing units to 
! ... store global variables, using PARALLEL_SUM MPI routine
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
      USE grid, ONLY: nij_l, myij
      USE pressure_epsilon, ONLY: p_g, pn_g, ep_g
      USE pressure_epsilon, ONLY: p, pn, ep
      USE th_capacity, ONLY: cp_g, ck_g
      USE th_capacity, ONLY: cp, ck
      USE gas_solid_viscosity, ONLY: mug_g, kapg_g
      USE gas_solid_viscosity, ONLY: mug, kapg
      USE turbulence, ONLY: scoeff,scoeff_g
!
      rgp_g   = 0.0D0
      rgpn_g = 0.0D0
      rog_g   = 0.0D0
      rlk_g   = 0.0D0
      rlkn_g  = 0.0D0
      p_g     = 0.0D0
      pn_g    = 0.0D0
      ep_g    = 0.0D0
      sieg_g  = 0.0D0
      siegn_g = 0.0D0
      tg_g    = 0.0D0
      siek_g  = 0.0D0
      siekn_g = 0.0D0
      ygc_g    = 0.0D0
      xgc_g    = 0.0D0
      rgpgc_g  = 0.0D0
      rgpgcn_g = 0.0D0
      tk_g    = 0.0D0
      ug_g    = 0.0D0
      vg_g    = 0.0D0
      uk_g    = 0.0D0
      vk_g    = 0.0D0
!
      cg_g    = 0.0
      cp_g    = 0.0
      ck_g    = 0.0
      mug_g   = 0.0
      kapg_g  = 0.0
      scoeff_g  = 0.0
!
      DO ij_l = 1, nij_l
        ij = myij(0,0,ij_l)
        rgp_g(ij) = rgp(ij_l)
        rgpn_g(ij) = rgpn(ij_l)
        rog_g(ij) = rog(ij_l)
        rlk_g(:,ij) = rlk(:,ij_l)
        rlkn_g(:,ij) = rlkn(:,ij_l)
        p_g(ij) = p(ij_l)
        pn_g(ij) = pn(ij_l)
        ep_g(ij) = ep(ij_l)
        sieg_g(ij) = sieg(ij_l)
        siegn_g(ij) = siegn(ij_l)
        tg_g(ij) = tg(ij_l)
        siek_g(:,ij) = siek(:,ij_l)
        siekn_g(:,ij) = siekn(:,ij_l)
        tk_g(:,ij) = tk(:,ij_l)
        rgpgcn_g(:, ij) = rgpgcn(:, ij_l)
        rgpgc_g(:, ij) = rgpgc(:, ij_l)
        ygc_g(:, ij) = ygc(:, ij_l)
        xgc_g(:, ij) = xgc(:, ij_l)
        ug_g(ij) = ug(ij_l)
        vg_g(ij) = vg(ij_l)
        uk_g(:,ij) = uk(:,ij_l)
        vk_g(:,ij) = vk(:,ij_l)
!
        ck_g(:,ij) = ck(:,ij_l)
        cp_g(:,ij) = cp(:,ij_l)
        cg_g(ij) = cg(ij_l)
        mug_g(ij) = mug(ij_l)
        kapg_g(ij) = kapg(ij_l)
        scoeff_g(ij) = scoeff(ij_l)
        
      END DO
!
      CALL parallel_sum_real(rgp_g, SIZE(rgp_g))
      CALL parallel_sum_real(rgpn_g, SIZE(rgpn_g))
      CALL parallel_sum_real(rog_g, SIZE(rog_g))
      CALL parallel_sum_real(rlk_g, SIZE(rlk_g))
      CALL parallel_sum_real(rlkn_g, SIZE(rlkn_g))
      CALL parallel_sum_real(p_g, SIZE(p_g))
      CALL parallel_sum_real(pn_g, SIZE(pn_g))
      CALL parallel_sum_real(ep_g, SIZE(ep_g) )
      CALL parallel_sum_real(sieg_g, SIZE(sieg_g) )
      CALL parallel_sum_real(siegn_g, SIZE(siegn_g) )
      CALL parallel_sum_real(tg_g, SIZE(tg_g) )
      CALL parallel_sum_real(siek_g, SIZE(siek_g) )
      CALL parallel_sum_real(siekn_g, SIZE(siekn_g) )
      CALL parallel_sum_real(ygc_g, SIZE(ygc_g))
      CALL parallel_sum_real(xgc_g, SIZE(xgc_g))
      CALL parallel_sum_real(rgpgcn_g, SIZE(rgpgcn_g))
      CALL parallel_sum_real(rgpgc_g, SIZE(rgpgc_g))
      CALL parallel_sum_real(tk_g, SIZE(tk_g) )
      CALL parallel_sum_real(ug_g, SIZE(ug_g) )
      CALL parallel_sum_real(vg_g, SIZE(vg_g) )
      CALL parallel_sum_real(uk_g, SIZE(uk_g) )
      CALL parallel_sum_real(vk_g, SIZE(vk_g) )
!
      CALL parallel_sum_real(cp_g, SIZE(cp_g) )
      CALL parallel_sum_real(ck_g, SIZE(ck_g) )
      CALL parallel_sum_real(cg_g, SIZE(cg_g) )
      CALL parallel_sum_real(mug_g, SIZE(mug_g) )
      CALL parallel_sum_real(kapg_g, SIZE(kapg_g) )
      CALL parallel_sum_real(scoeff_g, SIZE(scoeff_g) )
!
      END SUBROUTINE
