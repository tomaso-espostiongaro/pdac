!----------------------------------------------------------------------
      MODULE tilde_momentum
!----------------------------------------------------------------------
      IMPLICIT NONE
!
! ... tilde momentum density
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: rug, rvg, rwg
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rus, rvs, rws
      REAL*8, DIMENSION(:),   ALLOCATABLE :: rugn, rvgn, rwgn
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rusn, rvsn, rwsn
!
! ... convective momentum fluxes
!
      REAL*8, ALLOCATABLE :: ugfe(:), ugfn(:), ugft(:)
      REAL*8, ALLOCATABLE :: vgfe(:), vgfn(:), vgft(:)
      REAL*8, ALLOCATABLE :: wgfe(:), wgfn(:), wgft(:)
!
      REAL*8, ALLOCATABLE :: usfe(:,:), usfn(:,:), usft(:,:)
      REAL*8, ALLOCATABLE :: vsfe(:,:), vsfn(:,:), vsft(:,:)
      REAL*8, ALLOCATABLE :: wsfe(:,:), wsfn(:,:), wsft(:,:)
!
! ... gas-particle drag coefficient
!
      REAL*8, DIMENSION(:), ALLOCATABLE  :: kpgv
!
! ... momentum exchange coefficients
      REAL*8, DIMENSION(:,:), ALLOCATABLE  :: appu, appv, appw
      PUBLIC
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_momentum
      USE dimensions
      USE domain_mapping, ONLY: ncint, ncdom
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      IMPLICIT NONE
!
      ALLOCATE( rugn(ncint), rwgn(ncint))
      ALLOCATE( rusn(ncint,nsolid), rwsn(ncint,nsolid))
      rugn = 0.D0; rwgn = 0.D0
      rusn = 0.D0; rwsn = 0.D0
      IF (job_type == JOB_TYPE_3D) THEN
        ALLOCATE( rvgn(ncint) )
        ALLOCATE( rvsn(ncint,nsolid) )
        rvgn = 0.D0
        rvsn = 0.D0
      END IF
!
      RETURN
      END SUBROUTINE allocate_momentum
!-----------------------------------------------------
      SUBROUTINE allocate_fluxes
!
      USE dimensions, ONLY: nsolid
      USE domain_mapping, ONLY: ncint, ncdom
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE tilde_energy, ONLY: rhg, rhs
      IMPLICIT NONE
!
! ... Allocate and initialize "tilde" terms (gas)
!
      ALLOCATE(rug(ncdom), rwg(ncdom))
      rug = 0.0D0
      rwg = 0.0D0
      IF (job_type == JOB_TYPE_3D) THEN
        ALLOCATE(rvg(ncdom))
        rvg = 0.0D0
      END IF
      ALLOCATE(rhg(ncdom))
      rhg = 0.D0
!
! ... Allocate and initialize "tilde" terms (particles).
!
      ALLOCATE(rus(ncdom,nsolid), rws(ncdom,nsolid))
      rus = 0.0D0
      rws = 0.0D0
      IF (job_type == JOB_TYPE_3D) THEN
        ALLOCATE(rvs(ncdom,nsolid))
        rvs = 0.0D0
      END IF
      ALLOCATE(rhs(ncdom,nsolid))
      rhs = 0.D0
!
! ... Allocate and initialize interphase terms.
!
      ALLOCATE(appu(ncdom, ((nsolid+1)**2+(nsolid+1))/2),   &
               appw(ncdom, ((nsolid+1)**2+(nsolid+1))/2))

      appu = 0.0D0
      appw = 0.0D0
      IF (job_type == JOB_TYPE_3D) THEN
        ALLOCATE(appv(ncdom, ((nsolid+1)**2+(nsolid+1))/2) )
        appv = 0.0D0
      END IF
!
      RETURN
      END SUBROUTINE allocate_fluxes
!-----------------------------------------------------
      SUBROUTINE deallocate_fluxes

        USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
        USE tilde_energy, ONLY: rhg, rhs

        IMPLICIT NONE
        
        DEALLOCATE(rug, rwg)
        DEALLOCATE(rus, rws)
        DEALLOCATE(appu, appw)
        IF (job_type == JOB_TYPE_3D) THEN
          DEALLOCATE(rvg)
          DEALLOCATE(rvs)
          DEALLOCATE(appv)
        END IF
        DEALLOCATE(rhg)
        DEALLOCATE(rhs)
        
      END SUBROUTINE deallocate_fluxes
!-----------------------------------------------------
      SUBROUTINE fieldn
!
! ... Data_exchange +
! ... Compute explicitly and store all fields and physical parameters
! ... at time ndt for explicit time-advancement
! ... (2D-3D-Compliant)
!
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions
      USE domain_mapping, ONLY: ncint, myijk, ncdom
      USE domain_mapping, ONLY: meshinds, data_exchange
      USE eos_gas, ONLY: xgc, ygc
      USE gas_components, ONLY: rgpgcn
      USE gas_solid_density, ONLY: rgp, rgpn, rlk, rlkn, rog
      USE gas_solid_temperature, ONLY: sieg, siegn, sies, siesn, tg, ts
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_viscosity, ONLY: viscon, mug, kapg
      USE gas_solid_viscosity, ONLY: gas_viscosity, part_viscosity
      USE grid, ONLY: dz, dy, dx, flag
      USE grid, ONLY: indx, indy, indz
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE pressure_epsilon, ONLY: p, ep, pn, epn
      USE set_indexes, ONLY: first_subscr, ijke, ijkn, ijkt
      USE set_indexes, ONLY: ctu1_subscr, ctu2_subscr, ctu3_subscr
      USE flux_limiters, ONLY: ctu
      USE io_files, ONLY: testunit
!
      IMPLICIT NONE
!
      REAL*8 :: dxp, dyp, dzp, indxp, indyp, indzp
      REAL*8 :: rgp_e, rgp_n, rgp_t, rlk_e, rlk_n, rlk_t
      INTEGER :: i, j, k, ijk, imesh, is, ig, info
      LOGICAL :: compute
!
! ... Data_exchange of primary fields
!
      CALL data_exchange(ug)
      CALL data_exchange(wg)
      IF (job_type == JOB_TYPE_3D) CALL data_exchange(vg)

      CALL data_exchange(us)
      CALL data_exchange(ws)
      IF (job_type == JOB_TYPE_3D) CALL data_exchange(vs)

      CALL data_exchange(p)
      CALL data_exchange(ep)
      CALL data_exchange(rgp)
      CALL data_exchange(rlk)

      CALL data_exchange(sieg)
      CALL data_exchange(sies)
      CALL data_exchange(tg)
      CALL data_exchange(ts)
!
! ... Initialize momentum 'tilde' terms (explicit momentum fluxes)
!
      DO ijk = 1, ncint

       ! ... WARNING: in the old code this was done only in fluid cells,
       ! ... but it is now preferable to loop over the whole domain
       !
       !compute = BTEST(flag(ijk),0)
       !IF( compute ) THEN

          CALL meshinds(ijk,imesh,i,j,k)
          CALL first_subscr(ijk)
          IF (ctu > 0) CALL ctu1_subscr(ijk)
          IF (ctu > 1) CALL ctu2_subscr(ijk)
          IF (ctu > 2) CALL ctu3_subscr(ijk)
!
          IF (job_type == JOB_TYPE_2D ) THEN

            dxp = dx(i) + dx(i+1)
            dzp = dz(k) + dz(k+1)
            indxp = 1.D0 / dxp
            indzp = 1.D0 / dzp

            rgp_e = ( dx(i+1) * rgp(ijk) + dx(i) * rgp(ijke) ) * indxp
            rgp_t = ( dz(k+1) * rgp(ijk) + dz(k) * rgp(ijkt) ) * indzp
            rugn(ijk)  = rgp_e * ug(ijk)
            rwgn(ijk)  = rgp_t * wg(ijk)
           
            DO is = 1, nsolid
               rlk_e = ( rlk(ijk,is) * dx(i+1) + rlk(ijke,is) * dx(i) ) * indxp
               rlk_t = ( rlk(ijk,is) * dz(k+1) + rlk(ijkt,is) * dz(k) ) * indzp
               rusn(ijk,is)  = rlk_e * us(ijk,is)
               rwsn(ijk,is)  = rlk_t * ws(ijk,is)
            END DO
 
          ELSE IF (job_type == JOB_TYPE_3D) THEN

            dxp = dx(i) + dx(i+1)
            dyp = dy(j) + dy(j+1)
            dzp = dz(k) + dz(k+1)
            indxp = 1.D0 / dxp
            indyp = 1.D0 / dyp
            indzp = 1.D0 / dzp

            rgp_e = ( dx(i+1) * rgp(ijk) + dx(i) * rgp(ijke) ) * indxp
            rgp_n = ( dy(j+1) * rgp(ijk) + dy(j) * rgp(ijkn) ) * indyp
            rgp_t = ( dz(k+1) * rgp(ijk) + dz(k) * rgp(ijkt) ) * indzp
            rugn(ijk)  = rgp_e * ug(ijk)
            rvgn(ijk)  = rgp_n * vg(ijk)
            rwgn(ijk)  = rgp_t * wg(ijk)
!     
            DO is = 1, nsolid
              rlk_e = ( rlk(ijk,is) * dx(i+1) + rlk(ijke,is) * dx(i) ) * indxp
              rlk_n = ( rlk(ijk,is) * dy(j+1) + rlk(ijkn,is) * dy(j) ) * indyp
              rlk_t = ( rlk(ijk,is) * dz(k+1) + rlk(ijkt,is) * dz(k) ) * indzp
              rusn(ijk,is)  = rlk_e * us(ijk,is)
              rvsn(ijk,is)  = rlk_n * vs(ijk,is)
              rwsn(ijk,is)  = rlk_t * ws(ijk,is)
            END DO
         
          END IF 
!
! ... Store fields at time n*dt
!
          pn(ijk)    = p(ijk)
          epn(ijk)   = ep(ijk)
          rgpn(ijk)  = rgp(ijk)
          siegn(ijk) = sieg(ijk)
!
          DO ig = 1, ngas
            rgpgcn(ijk,ig) = rgpn(ijk) * ygc(ijk,ig)
          END DO
!
          DO is = 1, nsolid
            rlkn(ijk,is)  = rlk(ijk,is)
            siesn(ijk,is) = sies(ijk,is)
          END DO

        !END IF
        ! ...
!
! ... Compute the temperature-dependent gas viscosity and th. conductivity
! ... For inviscid simulation (gas_viscosity = FALSE) mug is used only
! ... to compute the gas-particle drag coefficient
!
         CALL viscon( mug(ijk), kapg(ijk), xgc(ijk,:), tg(ijk) )
      END DO
!
      CALL data_exchange(pn)
      CALL data_exchange(epn)
      CALL data_exchange(rgpn)
      CALL data_exchange(rlkn)

      RETURN
      END SUBROUTINE fieldn
!----------------------------------------------------------------------
      SUBROUTINE tilde
!
      USE atmospheric_conditions, ONLY: gravz, gravx, gravy
      USE dimensions
      USE domain_mapping, ONLY: meshinds
      USE domain_mapping, ONLY: ncint, myijk, ncdom, data_exchange
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE gas_solid_density, ONLY: rgp, rlk, rgpn, rlkn
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE grid, ONLY: dx, dy, dz, flag
      USE grid, ONLY: indx, indy, indz, inr, inrb
      USE momentum_transfer, ONLY: kdrags, inter, sdrag, drag_model
      USE pressure_epsilon, ONLY: ep, p
      USE time_parameters, ONLY: dt, time, alpha, alphagrav
      USE gas_solid_viscosity, ONLY: viscg, viscs
      USE gas_solid_viscosity, ONLY: gas_viscosity, part_viscosity
      USE gas_solid_viscosity, ONLY: mug
      USE gas_solid_viscosity, ONLY: gvisx, gvisy, gvisz, pvisx, pvisy, pvisz
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE immersed_boundaries, ONLY: fptx, fpty, fptz, numx, numy, numz
      USE immersed_boundaries, ONLY: immb
      USE mass_sink, ONLY: sink
      USE particles_constants, ONLY: inrl
      USE pressure_epsilon, ONLY: p, ep
      USE pressure_epsilon, ONLY: pn, epn, pmodel
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm, ijkt, ijke, ijkn
      USE set_indexes, ONLY: ctu1_subscr, ctu2_subscr, ctu3_subscr
      USE turbulence_model, ONLY: mugt
      USE flux_limiters, ONLY: ctu
!
      IMPLICIT NONE
!
      INTEGER :: i, j, k, is, imesh
      INTEGER :: ijk, ll, l, ls, ls1
      REAL*8 :: dxp, dyp, dzp, indxp, indyp, indzp
      REAL*8 :: dxi, dxip1, dyj, dyjp1, dzk, dzkp1
      REAL*8 :: ugfw, ugfs, ugfb, vgfw, vgfs, vgfb, wgfw, wgfs, wgfb
      REAL*8 :: ugfx, ugfy, ugfz, vgfx, vgfy, vgfz, wgfx, wgfy, wgfz
      REAL*8 :: usfw, usfs, usfb, vsfw, vsfs, vsfb, wsfw, wsfs, wsfb
      REAL*8 :: usfx, usfy, usfz, vsfx, vsfy, vsfz, wsfx, wsfy, wsfz
      REAL*8 :: rgp_e, rgp_n, rgp_t, rlk_e, rlk_n, rlk_t
      REAL*8 :: rug_tmp, rvg_tmp, rwg_tmp
      REAL*8 :: rus_tmp, rvs_tmp, rws_tmp
      REAL*8 :: epn_e, epn_n, epn_t, epsn_e, epsn_n, epsn_t
      REAL*8 :: force, presn, dragn
      REAL*8 :: pseudou, pseudov, pseudow, ep_e, ep_n, ep_t
      INTEGER :: fx, fy, fz
      REAL*8, ALLOCATABLE :: dugs(:), dvgs(:), dwgs(:)
      REAL*8, ALLOCATABLE :: nul(:)
      LOGICAL :: compute
!
! ... Allocate and initialize gas and particle viscous stresses
!
      ALLOCATE( gvisx(ncint), gvisz(ncint) )
      ALLOCATE( pvisx(ncint,nsolid), pvisz(ncint,nsolid) )
      gvisx = 0.D0; gvisz = 0.D0
      pvisx = 0.D0; pvisz = 0.D0

      IF ( job_type == JOB_TYPE_3D ) THEN
        ALLOCATE( gvisy(ncint) )
        ALLOCATE( pvisy(ncint,nsolid) )
        gvisy = 0.D0
        pvisy = 0.D0
      END IF
!
! ... Calculate gas viscous stress tensor
!
      IF ( gas_viscosity ) CALL viscg        
!
! ... Calculate particles viscous stress tensor
!
      IF ( part_viscosity ) CALL viscs
!
! ... Allocate and initialize gas convective fluxes
!
      ALLOCATE( ugfe(ncdom), ugft(ncdom) )
      ALLOCATE( wgfe(ncdom), wgft(ncdom) )
      ugfe = 0.0D0; ugft = 0.0D0
      ugfw = 0.0D0; ugfb = 0.0D0
      wgfe = 0.0D0; wgft = 0.0D0
      wgfw = 0.0D0; wgfb = 0.0D0

      IF (job_type == JOB_TYPE_3D) THEN
        ALLOCATE( ugfn(ncdom) )
        ALLOCATE( vgfe(ncdom), vgfn(ncdom), vgft(ncdom) )
        ALLOCATE( wgfn(ncdom) )
        ugfn = 0.0D0
        ugfs = 0.0D0
        vgfe = 0.0D0; vgfn = 0.0D0; vgft = 0.0D0
        vgfw = 0.0D0; vgfs = 0.0D0; vgfb = 0.0D0
        wgfn = 0.0D0
        wgfs = 0.0D0
      END IF
!
! ... Allocate and initialize particles convective fluxes
!
      ALLOCATE( usfe(ncdom,nsolid), usft(ncdom,nsolid) )
      ALLOCATE( wsfe(ncdom,nsolid), wsft(ncdom,nsolid) )

      usfe = 0.0D0; usft = 0.0D0
      usfw = 0.0D0; usfb = 0.0D0
      wsfe = 0.0D0; wsft = 0.0D0
      wsfw = 0.0D0; wsfb = 0.0D0

      IF (job_type == JOB_TYPE_3D) THEN
        ALLOCATE( usfn(ncdom,nsolid) )
        ALLOCATE( vsfe(ncdom,nsolid), vsfn(ncdom,nsolid), vsft(ncdom,nsolid) )
        ALLOCATE( wsfn(ncdom,nsolid) )
        usfn = 0.0D0
        usfs = 0.0D0
        vsfe = 0.0D0; vsfn = 0.0D0; vsft = 0.0D0
        vsfw = 0.0D0; vsfs = 0.0D0; vsfb = 0.0D0
        wsfn = 0.0D0
        wsfs = 0.0D0
      END IF
!
! ... Allocate and initialize gas-particle drag coefficient
!
      ALLOCATE( kpgv(nsolid) )
      ALLOCATE( dugs(nsolid) )
      ALLOCATE( dvgs(nsolid) )
      ALLOCATE( dwgs(nsolid) )
      kpgv = 0.0D0
      dugs = 0.0D0
      dvgs = 0.0D0
      dwgs = 0.0D0
!
! ... (a temporary array only used in 2D) ...
!
      ALLOCATE( nul( ( ( nsolid + 1 )**2 + ( nsolid + 1 ) ) / 2 ) )
      nul = 0.D0
!
! ... Compute all convective East, North, and Top fluxes 
! ... in the physical domain and ghost cells

      CALL compute_all_fluxes_tilde
!      CALL test_fluxes_tilde
!
! ... Fluxes on West, South and Bottom sides keep values 
! ... of East, North and Top fluxes from neighbouring cells.
!
      mesh_loop: DO ijk = 1, ncint
        compute = BTEST(flag(ijk),0)

        IF( compute ) THEN
          CALL meshinds(ijk,imesh,i,j,k)
          CALL subscr(ijk)
          IF (ctu > 0) CALL ctu1_subscr(ijk)
          IF (ctu > 1) CALL ctu2_subscr(ijk)
          IF (ctu > 2) CALL ctu3_subscr(ijk)
!
          ugfw = ugfe(imjk)
          ugfb = ugft(ijkm)
          wgfw = wgfe(imjk)
          wgfb = wgft(ijkm)

          ugfx = ugfe(ijk) - ugfw
          ugfz = ugft(ijk) - ugfb
          wgfx = wgfe(ijk) - wgfw
          wgfz = wgft(ijk) - wgfb
!
          IF (job_type == JOB_TYPE_2D ) THEN
            ugfy = 0.D0
            vgfx = 0.D0
            vgfy = 0.D0
            vgfz = 0.D0
            wgfy = 0.D0
          ELSE IF (job_type == JOB_TYPE_3D) THEN
            ugfs = ugfn(ijmk)
            vgfw = vgfe(imjk)
            vgfs = vgfn(ijmk)
            vgfb = vgft(ijkm)
            wgfs = wgfn(ijmk)

            ugfy = ugfn(ijk) - ugfs
            vgfx = vgfe(ijk) - vgfw
            vgfy = vgfn(ijk) - vgfs
            vgfz = vgft(ijk) - vgfb
            wgfy = wgfn(ijk) - wgfs
          END IF
!
! ... compute explicit (tilde) terms in the momentum equation (gas)
! 
          dxi   = dx(i)
          dxip1 = dx(i+1)
          dyj   = dy(j)
          dyjp1 = dy(j+1)
          dzk   = dz(k)
          dzkp1 = dz(k+1)
          dxp   = dx(i) + dx(i+1)
          dyp   = dy(j) + dy(j+1)
          dzp   = dz(k) + dz(k+1)
          indxp = 1.D0 / dxp
          indyp = 1.D0 / dyp
          indzp = 1.D0 / dzp

          IF ( pmodel == 1 ) THEN
            epn_e = (dxi*epn(ijke) + dxip1*epn(ijk)) * indxp
            epn_n = (dyj*epn(ijkn) + dyjp1*epn(ijk)) * indyp
            epn_t = (dzk*epn(ijkt) + dzkp1*epn(ijk)) * indzp
          ELSE IF ( pmodel == 2) THEN
            epn_e = 1.D0
            epn_n = 1.D0
            epn_t = 1.D0
          END IF
!         
          rgp_e = (dx(i+1)*rgp(ijk)+dx(i)*rgp(ijke)) * indxp
          rgp_n = (dy(j+1)*rgp(ijk)+dy(j)*rgp(ijkn)) * indyp
          rgp_t = (dz(k+1)*rgp(ijk)+dz(k)*rgp(ijkt)) * indzp
!
          rug_tmp = gvisx(ijk)                     
          rug_tmp = rug_tmp - indxp * 2.D0 * ugfx * inrb(i) 
          rug_tmp = rug_tmp + (1.D0-alphagrav)*(dx(i+1)*rgpn(ijk)+dx(i)*rgpn(ijke))*indxp * gravx  
          rug_tmp = rug_tmp - indy(j) * ugfy                  
          rug_tmp = rug_tmp - indz(k) * ugfz  
          rug_tmp = rug_tmp + (1.D0-alpha) * indxp *2.D0* epn_e * (pn(ijk)-pn(ijke)) 
          rug (ijk) = rugn(ijk) + dt * rug_tmp
          !
          IF (job_type == JOB_TYPE_3D) THEN
            rvg_tmp = gvisy(ijk)                     
            rvg_tmp = rvg_tmp + (1.D0-alphagrav)*(dy(j+1)*rgpn(ijk)+dy(j)*rgpn(ijkn))*indyp*gravy  
            rvg_tmp = rvg_tmp - indx(i) * vgfx               
            rvg_tmp = rvg_tmp - indyp * 2.D0 * vgfy   
            rvg_tmp = rvg_tmp - indz(k) * vgfz    
            rvg_tmp = rvg_tmp + (1.D0-alpha) * indyp *2.D0* epn_n * (pn(ijk)-pn(ijkn)) 
            rvg(ijk) = rvgn(ijk) + dt * rvg_tmp
            !
          END IF
!
          rwg_tmp = gvisz(ijk)                     
          rwg_tmp = rwg_tmp + (1.D0-alphagrav)*(dz(k+1)*rgpn(ijk)+dz(k)*rgpn(ijkt))*indzp * gravz  
          rwg_tmp = rwg_tmp - indx(i) * wgfx * inr(i)
          rwg_tmp = rwg_tmp - indy(j) * wgfy
          rwg_tmp = rwg_tmp - indzp * 2.D0 * wgfz
          rwg_tmp = rwg_tmp + (1.D0-alpha) * indzp *2.D0* epn_t * (pn(ijk)-pn(ijkt)) 
          rwg(ijk) = rwgn(ijk) + dt * rwg_tmp
          !
!
! ... same procedure carried out for each particulate phases
!
          DO is = 1, nsolid
!
! ... West, South and Bottom fluxes (particles)
!
            usfw = usfe(imjk,is)
            usfb = usft(ijkm,is)
            wsfw = wsfe(imjk,is)
            wsfb = wsft(ijkm,is)

            usfx = usfe(ijk,is) - usfw
            usfz = usft(ijk,is) - usfb
            wsfx = wsfe(ijk,is) - wsfw
            wsfz = wsft(ijk,is) - wsfb

            IF (job_type == JOB_TYPE_2D ) THEN
              usfy = 0.D0
              vsfx = 0.D0
              vsfy = 0.D0
              vsfz = 0.D0
              wsfy = 0.D0
            ELSE IF (job_type == JOB_TYPE_3D) THEN
              usfs = usfn(ijmk,is)
              vsfw = vsfe(imjk,is)
              vsfs = vsfn(ijmk,is)
              vsfb = vsft(ijkm,is)
              wsfs = wsfn(ijmk,is)             

              usfy = usfn(ijk,is) - usfs
              vsfx = vsfe(ijk,is) - vsfw
              vsfy = vsfn(ijk,is) - vsfs
              vsfz = vsft(ijk,is) - vsfb
              wsfy = wsfn(ijk,is) - wsfs
            END IF

            IF ( pmodel == 1) THEN
              epsn_e = (dxi*rlkn(ijke,is) + dxip1*rlkn(ijk,is)) * indxp * inrl(is)
              epsn_n = (dyj*rlkn(ijkn,is) + dyjp1*rlkn(ijk,is)) * indyp * inrl(is)
              epsn_t = (dzk*rlkn(ijkt,is) + dzkp1*rlkn(ijk,is)) * indzp * inrl(is)
            ELSE IF ( pmodel == 2 ) THEN
              epsn_e = 0.D0
              epsn_n = 0.D0
              epsn_t = 0.D0
            END IF
!
! ... compute explicit (tilde) terms in the momentum equation (particles)
! 
            rlk_e = (dx(i+1)*rlk(ijk,is)+dx(i)*rlk(ijke,is)) * indxp
            rlk_n = (dy(j+1)*rlk(ijk,is)+dy(j)*rlk(ijkn,is)) * indyp
            rlk_t = (dz(k+1)*rlk(ijk,is)+dz(k)*rlk(ijkt,is)) * indzp
!
            rus_tmp = pvisx(ijk,is)               
            rus_tmp = rus_tmp - indxp * 2.D0 * usfx * inrb(i)
            rus_tmp = rus_tmp + (1.D0-alphagrav) * indxp * gravx * &
                      (dx(i+1)*rlkn(ijk,is)+dx(i)*rlkn(ijke,is))
            rus_tmp = rus_tmp - indy(j) * usfy  
            rus_tmp = rus_tmp - indz(k) * usfz 
            rus_tmp = rus_tmp + (1.D0-alpha) * indxp *2.D0* epsn_e * (pn(ijk)-pn(ijke)) 
            rus(ijk,is) = rusn(ijk,is) + dt * rus_tmp
            !
            ! ... add sink term (at bottom boundaries)
            rus(ijk,is) = rus(ijk,is) + dt * sink(ijk,is)*us(ijk,is)
!
            IF (job_type == JOB_TYPE_3D) THEN
              rvs_tmp = pvisy(ijk,is)              
              rvs_tmp = rvs_tmp + (1.D0-alphagrav) * indyp * gravy * &
                        (dy(j+1)*rlkn(ijk,is)+dy(j)*rlkn(ijkn,is))
              rvs_tmp = rvs_tmp - indx(i) * vsfx       
              rvs_tmp = rvs_tmp - indyp * 2.D0 * vsfy              
              rvs_tmp = rvs_tmp - indz(k) * vsfz    
              rvs_tmp = rvs_tmp + (1.D0-alpha) * indyp *2.D0* epsn_n * (pn(ijk)-pn(ijkn)) 
              rvs(ijk,is) = rvsn(ijk,is) + dt * rvs_tmp
              !
              ! ... add sink term (at bottom boundaries)
              rvs(ijk,is) = rvs(ijk,is) + dt * sink(ijk,is) * vs(ijk,is)
              !
            END IF
!
            rws_tmp = pvisz(ijk,is) 
            rws_tmp = rws_tmp + (1.D0-alphagrav) * indzp * gravz * &
                      ( rlkn(ijk,is) * dz(k+1) + rlkn(ijkt,is) * dz(k) )
            rws_tmp = rws_tmp - indx(i) * wsfx * inr(i)                 
            rws_tmp = rws_tmp - indy(j) * wsfy                      
            rws_tmp = rws_tmp - indzp * 2.D0 * wsfz
            rws_tmp = rws_tmp + (1.D0-alpha) * indzp *2.D0* epsn_t * (pn(ijk)-pn(ijkt)) 
            rws(ijk,is) = rwsn(ijk,is) + dt * rws_tmp
            !
            ! ... add sink term (at bottom boundaries)
            rws(ijk,is) = rws(ijk,is) + dt * sink(ijk,is) * ws(ijk,is)
            !
!
! ... Compute the gas-particle drag coefficients in the cell center
!
            dugs(is) = ( (ug(ijk)-us(ijk,is)) + (ug(imjk)-us(imjk,is)) )*0.5D0
            dwgs(is) = ( (wg(ijk)-ws(ijk,is)) + (wg(ijkm)-ws(ijkm,is)) )*0.5D0
            IF (job_type == JOB_TYPE_2D ) THEN
              dvgs(is) = 0.D0
            ELSE IF (job_type == JOB_TYPE_3D) THEN
              dvgs(is) = ( (vg(ijk)-vs(ijk,is)) + (vg(ijmk)-vs(ijmk,is)) )*0.5D0
            END IF
!
            IF (drag_model == 1) THEN
              CALL kdrags(kpgv(is), dugs(is), dvgs(is), dwgs(is), ep(ijk),     &
                    rgp(ijk), rlk(ijk,is), mug(ijk), is)                  
            ELSE IF (drag_model == 2) THEN
              CALL sdrag(kpgv(is), dugs(is), dvgs(is), dwgs(is), ep(ijk),     &
                    rgp(ijk), rlk(ijk,is), mug(ijk), is)                  
            END IF
          END DO
!
!
! ... Compute the particle-particle coefficients and the interphase matrix
!
          IF (job_type == JOB_TYPE_2D ) THEN
            CALL inter(appu(ijk,:), nul(:), appw(ijk,:), kpgv(:),    &
     &                 us, us, ws, rlk, ijk)
          ELSE IF (job_type == JOB_TYPE_3D) THEN
            CALL inter(appu(ijk,:), appv(ijk,:), appw(ijk,:), kpgv(:),    &
     &                 us, vs, ws, rlk, ijk)
          END IF

          rug(ijk) = rug(ijk) - (1.D0-alpha)*(dxi*appu(ijke,1)+dxip1*appu(ijk,1)) &
                                         *indxp*ug(ijk)
          IF (job_type == JOB_TYPE_3D) THEN
            rvg(ijk) = rvg(ijk) - (1.D0-alpha)*(dyj*appv(ijkn,1)+dyjp1*appv(ijk,1)) &
                                           *indyp*vg(ijk)
          END IF
          rwg(ijk) = rwg(ijk) - (1.D0-alpha)*(dzk*appw(ijkt,1)+dzkp1*appw(ijk,1)) &
                                         *indzp*wg(ijk)

      DO l=2,nphase

        ls1=l*(l-1)/2

        DO ll=1,l

          ls=ls1+ll

          IF (ll == 1) THEN

            rug(ijk) = rug(ijk) - (1.D0-alpha)*(dxi*appu(ijke,ls)+dxip1*appu(ijk,ls)) &
                                           *indxp*us(ijk,l-1)
            IF (job_type == JOB_TYPE_3D) THEN
              rvg(ijk) = rvg(ijk) - (1.D0-alpha)*(dyj*appv(ijkn,ls)+dyjp1*appv(ijk,ls)) & 
                                             *indyp*vs(ijk,l-1)
            END IF
            rwg(ijk) = rwg(ijk) - (1.D0-alpha)*(dzk*appw(ijkt,ls)+dzkp1*appw(ijk,ls)) &
                                           *indzp*ws(ijk,l-1)

            rus(ijk,l-1) = rus(ijk,l-1) - (1.D0-alpha)*(dxi*appu(ijke,ls)+dxip1*appu(ijk,ls)) &
                                                   *indxp*ug(ijk)
            IF (job_type == JOB_TYPE_3D) THEN
              rvs(ijk,l-1) = rvs(ijk,l-1) - (1.D0-alpha)*(dyj*appv(ijkn,ls)+dyjp1*appv(ijk,ls)) &
                                                     *indyp*vg(ijk)
            END IF
            rws(ijk,l-1) = rws(ijk,l-1) - (1.D0-alpha)*(dzk*appw(ijkt,ls)+dzkp1*appw(ijk,ls)) &
                                                   *indzp*wg(ijk)

          ELSE IF (ll == l) THEN

            rus(ijk,l-1) = rus(ijk,l-1) - (1.D0-alpha)*(dxi*appu(ijke,ls)+dxip1*appu(ijk,ls)) &
                                                   *indxp*us(ijk,l-1)
            IF (job_type == JOB_TYPE_3D) THEN
              rvs(ijk,l-1) = rvs(ijk,l-1) - (1.D0-alpha)*(dyj*appv(ijkn,ls)+dyjp1*appv(ijk,ls)) &
                                                     *indyp*vs(ijk,l-1)
            END IF
            rws(ijk,l-1) = rws(ijk,l-1) - (1.D0-alpha)*(dzk*appw(ijkt,ls)+dzkp1*appw(ijk,ls)) &
                                                  *indzp*ws(ijk,l-1)

          ELSE

            rus(ijk,l-1)  = rus(ijk,l-1)  - (1.D0-alpha)*(dxi*appu(ijke,ls)+dxip1*appu(ijk,ls)) &
                                                     *indxp*us(ijk,ll-1)
            rus(ijk,ll-1) = rus(ijk,ll-1) - (1.D0-alpha)*(dxi*appu(ijke,ls)+dxip1*appu(ijk,ls)) &
                                                     *indxp*us(ijk,l-1)
            IF (job_type == JOB_TYPE_3D) THEN
              rvs(ijk,l-1)  = rvs(ijk,l-1)  - (1.D0-alpha)*(dyj*appv(ijkn,ls)+dyjp1*appv(ijk,ls)) &
                                                       *indyp*vs(ijk,ll-1)
              rvs(ijk,ll-1) = rvs(ijk,ll-1) - (1.D0-alpha)*(dyj*appv(ijkn,ls)+dyjp1*appv(ijk,ls)) &
                                                       *indyp*vs(ijk,l-1)
            END IF
            rws(ijk,l-1)  = rws(ijk,l-1)  - (1.D0-alpha)*(dzk*appw(ijkt,ls)+dzkp1*appw(ijk,ls)) &
                                                     *indzp*ws(ijk,ll-1)
            rws(ijk,ll-1) = rws(ijk,ll-1) - (1.D0-alpha)*(dzk*appw(ijkt,ls)+dzkp1*appw(ijk,ls)) &
                                                     *indzp*ws(ijk,l-1)

          END IF

        END DO

      END DO
!
! ... On the immersed boundary the explicit terms must be modified
! ... by adding a force to mimic the boundary
!
          IF (immb == 1) THEN
            
            fx = numx(ijk)

            IF (fx/=0) THEN
              pseudou = fptx(fx)%vel(1)

              force = ( indxp * (dx(i+1)*rgp(ijk)+dx(i)*rgp(ijke)) * &
                        pseudou - rug(ijk) ) / dt
              ep_e  = (dx(i)*ep(ijke) + dx(i+1)*ep(ijk)) * indxp
              presn = - indxp * 2.D0 * ep_e * (p(ijke)-p(ijk))
              dragn = 0.D0
              DO is = 1, nsolid
                dragn = dragn - kpgv(is) * dugs(is)
              END DO
              force = force - presn - dragn
              rug(ijk) = rug(ijk) + dt * force

              DO is = 1, nsolid
                pseudou = fptx(fx)%vel(1+is)
                force = ( indxp * (dx(i+1)*rlk(ijk,is)+dx(i)*rlk(ijke,is)) * &
                          pseudou - rus(ijk,is) ) / dt
                ep_e  = (dx(i)*rlk(ijke,is) + dx(i+1)*rlk(ijk,is)) * indxp * &
                        inrl(is)
                presn = - indxp * 2.D0 * ep_e * (p(ijke)-p(ijk))
                dragn = kpgv(is) * dugs(is)
                force = force - presn - dragn
                rus(ijk,is) = rus(ijk,is) + dt * force
              END DO
            END IF
            
            IF (job_type == JOB_TYPE_3D) THEN

              fy = numy(ijk)

              IF (fy/=0) THEN
                pseudov = fpty(fy)%vel(1)
  
                force = ( indyp * (dy(j+1)*rgp(ijk)+dy(j)*rgp(ijkn)) * &
                          pseudov - rvg(ijk) ) / dt
                ep_n  = (dy(j)*ep(ijkn) + dy(j+1)*ep(ijk)) * indyp
                presn = - indyp * 2.D0 * ep_n * (p(ijkn)-p(ijk))
                dragn = 0.D0
                DO is = 1, nsolid
                  dragn = dragn - kpgv(is) * dvgs(is)
                END DO
                force = force - presn - dragn
                rvg(ijk) = rvg(ijk) + dt * force
 
                DO is = 1, nsolid
                  pseudov = fpty(fy)%vel(1+is)
                  force =( indyp * (dy(j+1)*rlk(ijk,is)+dy(j)*rlk(ijkn,is)) * &
                            pseudov - rvs(ijk,is) ) / dt
                  ep_n  =(dy(j)*rlk(ijkn,is) + dy(j+1)*rlk(ijk,is)) * indyp * &
                          inrl(is)
                  presn = - indyp * 2.D0 * ep_n * (p(ijkn)-p(ijk))
                  dragn = kpgv(is) * dvgs(is)
                  force = force - presn - dragn
                  rvs(ijk,is) = rvs(ijk,is) + dt * force
                END DO
              END IF
            END IF
            
            fz = numz(ijk)

            IF (fz/=0) THEN
              pseudow = fptz(fz)%vel(1)

              force = ( indzp * (dz(k+1)*rgp(ijk)+dz(k)*rgp(ijkt)) * &
                        pseudow - rwg(ijk) ) / dt
              ep_t = (dz(k)*ep(ijkt) + dz(k+1)*ep(ijk)) * indzp
              presn = - indzp * 2.D0 * ep_t * (p(ijkt)-p(ijk))
              dragn = 0.D0
              DO is = 1, nsolid
                dragn = dragn - kpgv(is) * dwgs(is)
              END DO
              force = force - presn - dragn
              rwg(ijk) = rwg(ijk) + dt * force

              DO is = 1, nsolid
                pseudow = fptz(fz)%vel(1+is)
                force = ( indzp * (dz(k+1)*rlk(ijk,is)+dz(k)*rlk(ijkt,is)) * &
                          pseudow - rws(ijk,is) ) / dt
                ep_t = (dz(k)*rlk(ijkt,is) + dz(k+1)*rlk(ijk,is)) * indzp * &
                       inrl(is) 
                presn = - indzp * 2.D0 * ep_t * (p(ijkt)-p(ijk))
                dragn = kpgv(is) * dwgs(is)
                force = force - presn - dragn
                rws(ijk,is) = rws(ijk,is) + dt * force
              END DO
            END IF

          END IF
          
        END IF
      END DO mesh_loop
!
      DEALLOCATE(ugfe, ugft)
      DEALLOCATE(wgfe, wgft)
      DEALLOCATE(gvisx, gvisz)
!
      DEALLOCATE(usfe, usft)
      DEALLOCATE(wsfe, wsft)
      DEALLOCATE(pvisx, pvisz)
!
      DEALLOCATE(kpgv)
      DEALLOCATE( dugs )
      DEALLOCATE( dvgs )
      DEALLOCATE( dwgs )
      DEALLOCATE(nul)
!
      CALL data_exchange(appu)
      CALL data_exchange(appw)
      CALL data_exchange(rug)
      CALL data_exchange(rwg)
      CALL data_exchange(rus)
      CALL data_exchange(rws)
!
      IF (job_type == JOB_TYPE_3D) THEN

        DEALLOCATE( ugfn )
        DEALLOCATE( vgfe, vgfn, vgft)
        DEALLOCATE( wgfn  )
        DEALLOCATE( gvisy )

        DEALLOCATE( usfn )
        DEALLOCATE( vsfe, vsfn, vsft)
        DEALLOCATE( wsfn  )
        DEALLOCATE( pvisy )

        CALL data_exchange(appv)
        CALL data_exchange(rvg)
        CALL data_exchange(rvs)

      END IF

      RETURN
      END SUBROUTINE tilde
!----------------------------------------------------------------------
      SUBROUTINE compute_all_fluxes_tilde
!
      USE dimensions, ONLY: nsolid, nx, ny, nz
      USE domain_mapping, ONLY: meshinds
      USE domain_mapping, ONLY: ncint, myijk, data_exchange
      USE convective_fluxes_u, ONLY: flu, muscl_flu, ctu1_flu, ctu2_flu, ctu3_flu
      USE convective_fluxes_v, ONLY: flv, muscl_flv, ctu1_flv, ctu2_flv, ctu3_flv
      USE convective_fluxes_w, ONLY: flw, muscl_flw, ctu1_flw, ctu2_flw, ctu3_flw
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE flux_limiters, ONLY: muscl, ctu
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE grid, ONLY: flag, fluid, bl_cell, immb_cell
      USE immersed_boundaries, ONLY: immb, b_e_, b_n_, b_t_
      USE interpolate_fields, ONLY: interpolate_x, interpolate_y, interpolate_z
      USE pressure_epsilon, ONLY: ep, p
      USE set_indexes, ONLY: subscr, stencil, ctu1_subscr, ctu2_subscr, ctu3_subscr
      USE set_indexes, ONLY: maskval, maskstencil, masks
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes, ONLY: nb, rnb
      USE set_indexes, ONLY: ctu1_nb, ctu1_rnb, ctu2_nb, ctu2_rnb, ctu3_nb, ctu3_rnb
!
      IMPLICIT NONE
!
      INTEGER :: ijk
      INTEGER :: i, j, k, is, imesh
      TYPE(stencil) :: u, v, w, dens
      TYPE(masks) :: mask_x, mask_y, mask_z
      TYPE(stencil) :: dens_stagx, dens_stagy, dens_stagz
      LOGICAL :: compute
!
! ... Compute fluxes on East, North and Top sides of a cell
! ... in the whole computational domain.
!
      DO ijk = 1, ncint
        compute  = BTEST(flag(ijk),0)
        !
        IF( compute ) THEN
          CALL subscr(ijk)
          IF (ctu > 0) CALL ctu1_subscr(ijk)
          IF (ctu > 1) CALL ctu2_subscr(ijk)
          IF (ctu > 2) CALL ctu3_subscr(ijk)
          CALL meshinds(ijk,imesh,i,j,k)
!
          IF (job_type == JOB_TYPE_2D ) THEN
            !
            ! ... Stencil Mask for Immersed Boundaries
            !
            IF (immb >= 1) THEN
              mask_x = maskval(1)
              mask_z = maskval(1)
              IF (flag(ijk)==bl_cell .OR. flag(ijk)==immb_cell) THEN
                !
                ! ... Stencil masks
                CALL rnb(mask_x,b_e_,ijk)
                CALL rnb(mask_z,b_t_,ijk)
              END IF
            END IF
!
! ... (GAS) ...
!
! ... Compute convective fluxes by using First Order Upwind
!
            ! ... Assemble computational stencil
            CALL rnb(u,ug,ijk)
            CALL rnb(w,wg,ijk)
            CALL nb(dens,rgp,ijk)
            IF (ctu > 0) THEN
              CALL ctu1_rnb(u,ug,ijk)
              CALL ctu1_rnb(w,wg,ijk)
              CALL ctu1_nb(dens,rgp,ijk)
              IF (ctu > 1) THEN
                CALL ctu2_rnb(u,ug,ijk)
                CALL ctu2_rnb(w,wg,ijk)
                CALL ctu2_nb(dens,rgp,ijk)
                IF (ctu > 2) THEN
                  CALL ctu3_rnb(u,ug,ijk)
                  CALL ctu3_rnb(w,wg,ijk)
                  CALL ctu3_nb(dens,rgp,ijk)
                END IF
              END IF
            END IF
!
            ! ... Mask non-physical velocities from Immersed Boundaries
            !
            IF (immb >= 1) THEN
              IF (flag(ijk)==bl_cell .OR. flag(ijk)==immb_cell) THEN
                u = maskstencil(u,mask_x)
                w = maskstencil(w,mask_z)
              END IF
            END IF
!
            ! ... Interpolate density on the staggered grid
            CALL interpolate_x(dens, dens_stagx, i)
            CALL interpolate_z(dens, dens_stagz, k)

            CALL flu(ugfe(ijk), ugft(ijk), ugfe(imjk), ugft(ijkm),  & 
                     dens_stagx, u, w, i)

            CALL flw(wgfe(ijk), wgft(ijk), wgfe(imjk), wgft(ijkm),  &
                     dens_stagz, u, w, i, k)
!
! ... Second order MUSCL correction
!
            IF (muscl > 0 .AND. flag(ijk)==fluid) THEN
              IF ( i /= nx-1 ) CALL muscl_flu(ugfe(ijk), ugft(ijk),   &
                                        dens_stagx, u, w, i, k)
              IF ( k /= nz-1 ) CALL muscl_flw(wgfe(ijk), wgft(ijk),   &
                                        dens_stagz, u, w, i, k)
            END IF

            IF (ctu > 0) THEN
!
! ... First order Corner Transport Upwind correction (step 2)
!
              CALL ctu1_flu(ugfe(ijk), ugft(ijk), dens_stagx, u, w, i, k)
              CALL ctu1_flw(wgfe(ijk), wgft(ijk), dens_stagz, u, w, i, k)
!
! ... Second order Corner Transport Upwind correction (step 3)
!
              IF (ctu > 1 .AND. muscl == 0 .AND. flag(ijk) == fluid) THEN
                IF (i /= nx-1) CALL ctu2_flu(ugfe(ijk), ugft(ijk), &
                                         dens_stagx, u, w, i, k)
                IF (k /= nz-1) CALL ctu2_flw(wgfe(ijk), wgft(ijk), &
                                         dens_stagz, u, w, i, k)
!
! ... Third order Corner Transport Upwind correction (step 4)
!
                IF (ctu > 2) THEN
                  CALL ctu3_flu(ugfe(ijk), ugft(ijk), dens_stagx, u, w, i, k)
                  CALL ctu3_flw(wgfe(ijk), wgft(ijk), dens_stagz, u, w, i, k)
                END IF

              END IF
            END IF

! ... (PARTICLES) ...
!
! ... Compute convective fluxes by using First Order Upwind
!
            DO is = 1, nsolid

              ! ... Assemble computational stencil
              CALL rnb(u,us(:,is),ijk)
              CALL rnb(w,ws(:,is),ijk)
              CALL nb(dens,rlk(:,is),ijk)

              IF (ctu > 0) THEN
                CALL ctu1_rnb(u,us(:,is),ijk)
                CALL ctu1_rnb(w,ws(:,is),ijk)
                CALL ctu1_nb(dens,rlk(:,is),ijk)
                IF (ctu > 1) THEN
                  CALL ctu2_rnb(u,us(:,is),ijk)
                  CALL ctu2_rnb(w,ws(:,is),ijk)
                  CALL ctu2_nb(dens,rlk(:,is),ijk)
                  IF (ctu > 2) THEN
                    CALL ctu3_rnb(u,us(:,is),ijk)
                    CALL ctu3_rnb(w,ws(:,is),ijk)
                    CALL ctu3_nb(dens,rlk(:,is),ijk)
                  END IF
                END IF
              END IF

              ! ... Mask non-physical velocities from Immersed Boundaries
              !
              IF (immb >= 1) THEN
                IF (flag(ijk)==bl_cell .OR. flag(ijk)==immb_cell) THEN
                  u = maskstencil(u,mask_x)
                  w = maskstencil(w,mask_z)
                END IF
              END IF

              ! ... Interpolate density on the staggered grid
              CALL interpolate_x(dens, dens_stagx, i)
              CALL interpolate_z(dens, dens_stagz, k)

              CALL flu(usfe(ijk,is), usft(ijk,is), usfe(imjk,is), usft(ijkm,is), &
                       dens_stagx, u, w, i)

              CALL flw(wsfe(ijk,is), wsft(ijk,is), wsfe(imjk,is), wsft(ijkm,is), &
                       dens_stagz, u, w, i, k)
!
! ... Second order MUSCL correction
!
              IF (muscl > 0 .AND. flag(ijk)==fluid) THEN
                IF ( i /= nx-1 ) CALL muscl_flu(usfe(ijk,is), usft(ijk,is),  &
                                          dens_stagx, u, w, i, k)
                IF ( k /= nz-1 ) CALL muscl_flw(wsfe(ijk,is), wsft(ijk,is),  &
                                          dens_stagz, u, w, i, k)
              END IF

              IF (ctu > 0) THEN
!
! ... First order Corner Transport Upwind correction (step 2)
!
                CALL ctu1_flu(usfe(ijk,is), usft(ijk,is), dens_stagx, u, w, i, k)
                CALL ctu1_flw(wsfe(ijk,is), wsft(ijk,is), dens_stagz, u, w, i, k)
!
! ... Second order Corner Transport Upwind correction (step 3)
!
                IF (ctu > 1 .AND. muscl == 0 .AND. flag(ijk) == fluid) THEN
                  IF (i/=nx-1) CALL ctu2_flu(usfe(ijk,is), usft(ijk,is), &
                                      dens_stagx, u, w, i, k)
                  IF (k/=nz-1) CALL ctu2_flw(wsfe(ijk,is), wsft(ijk,is), &
                                           dens_stagz, u, w, i, k)
!
! ... Third order Corner Transport Upwind correction (step4)
!
                  IF (ctu > 2) THEN
                    CALL ctu3_flu(usfe(ijk,is), usft(ijk,is), dens_stagx, u, w, i, k)
                    CALL ctu3_flw(wsfe(ijk,is), wsft(ijk,is), dens_stagz, u, w, i, k)
                  END IF
                END IF
              END IF

            END DO
            
          ELSE IF (job_type == JOB_TYPE_3D) THEN
            !
            ! ... Stencil Mask for Immersed Boundaries
            !
            IF (immb >= 1) THEN
              mask_x = maskval(1)
              mask_y = maskval(1)
              mask_z = maskval(1)
              IF (flag(ijk)==bl_cell .OR. flag(ijk)==immb_cell) THEN
                !
                ! ... Stencil masks
                CALL rnb(mask_x,b_e_,ijk)
                CALL rnb(mask_y,b_n_,ijk)
                CALL rnb(mask_z,b_t_,ijk)
              END IF
            END IF
!
! ... (GAS) ...
!
! ... Compute convective fluxes by using First Order Upwind
!
            ! ... Assemble computational stencil
            CALL rnb(u,ug,ijk)
            CALL rnb(v,vg,ijk)
            CALL rnb(w,wg,ijk) 
            CALL nb ( dens, rgp, ijk )
!
            IF (ctu > 0) THEN
              CALL ctu1_rnb(u,ug,ijk)
              CALL ctu1_rnb(v,vg,ijk)
              CALL ctu1_rnb(w,wg,ijk)
              CALL ctu1_nb(dens,rgp,ijk)
              IF (ctu > 1) THEN
                CALL ctu2_rnb(u,ug,ijk)
                CALL ctu2_rnb(v,vg,ijk)
                CALL ctu2_rnb(w,wg,ijk)
                CALL ctu2_nb(dens,rgp,ijk)
                IF (ctu > 2) THEN
                  CALL ctu3_rnb(u,ug,ijk)
                  CALL ctu3_rnb(v,vg,ijk)
                  CALL ctu3_rnb(w,wg,ijk)
                  CALL ctu3_nb(dens,rgp,ijk)
                END IF
              END IF
            END IF
!

            ! ... Mask non-physical velocities from Immersed Boundaries
            !
            IF (immb >= 1) THEN
              IF (flag(ijk)==bl_cell .OR. flag(ijk)==immb_cell) THEN
                u = maskstencil(u,mask_x)
                v = maskstencil(v,mask_y)
                w = maskstencil(w,mask_z)
              END IF
            END IF

            ! ... Interpolate density on the staggered grid
            CALL interpolate_x(dens, dens_stagx, i)
            CALL interpolate_y(dens, dens_stagy, j)
            CALL interpolate_z(dens, dens_stagz, k)

            CALL flu(ugfe(ijk), ugfn(ijk), ugft(ijk),                        &
                     ugfe(imjk), ugfn(ijmk), ugft(ijkm),                     &
                     dens_stagx, u, v, w, i)

            CALL flv(vgfe(ijk), vgfn(ijk), vgft(ijk),                        &
                     vgfe(imjk), vgfn(ijmk), vgft(ijkm),                     &
                     dens_stagy, u, v, w, j)

            CALL flw(wgfe(ijk), wgfn(ijk), wgft(ijk),                        &
                     wgfe(imjk), wgfn(ijmk), wgft(ijkm),                     &
                     dens_stagz, u, v, w, k)
!
! ... Second order MUSCL correction
!
            IF (muscl > 0 .AND. flag(ijk)==fluid) THEN
              IF ( i /= nx-1 ) &
                CALL muscl_flu(ugfe(ijk), ugfn(ijk), ugft(ijk), &
                         dens_stagx, u, v, w, i, j, k)
              IF ( j /= ny-1 ) &
                CALL muscl_flv(vgfe(ijk), vgfn(ijk), vgft(ijk), &
                         dens_stagy, u, v, w, i, j, k)
              IF ( k /= nz-1 ) &
                CALL muscl_flw(wgfe(ijk), wgfn(ijk), wgft(ijk), &
                         dens_stagz, u, v, w, i, j, k)
            END IF
!
            IF (ctu > 0) THEN
!
! ... First order Corner Transport Upwind correction
!
              CALL ctu1_flu(ugfe(ijk), ugfn(ijk), ugft(ijk), &
                        dens_stagx, u, v, w, i, j, k)
              CALL ctu1_flv(vgfe(ijk), vgfn(ijk), vgft(ijk), &
                        dens_stagy, u, v, w, i, j, k)
              CALL ctu1_flw(wgfe(ijk), wgfn(ijk), wgft(ijk), &
                        dens_stagz, u, v, w, i, j, k)
!
! ... Second order Corner Transport Upwind correction
!
              IF (ctu > 1 .AND. muscl == 0 .AND. flag(ijk) == fluid) THEN
                IF (i /= nx-1) CALL ctu2_flu(ugfe(ijk), ugfn(ijk), ugft(ijk), &
                                         dens_stagx, u, v, w, i, j, k)
                IF (j /= ny-1) CALL ctu2_flv(vgfe(ijk), vgfn(ijk), vgft(ijk), &
                                         dens_stagy, u, v, w, i, j, k)
                IF (k /= nz-1) CALL ctu2_flw(wgfe(ijk), wgfn(ijk), wgft(ijk), &
                                         dens_stagz, u, v, w, i, j, k)
!
! ... Third order Corner Transport Upwind correction
!
                IF (ctu > 2) THEN
                  CALL ctu3_flu(ugfe(ijk), ugfn(ijk), ugft(ijk), dens_stagx, u, v, w, i, j, k)
                  CALL ctu3_flv(vgfe(ijk), vgfn(ijk), vgft(ijk), dens_stagy, u, v, w, i, j, k)
                  CALL ctu3_flw(wgfe(ijk), wgfn(ijk), wgft(ijk), dens_stagz, u, v, w, i, j, k)
                END IF
              END IF
            END IF

!
! ... (PARTICLES) ...
!
! ... Compute convective fluxes by using First Order Upwind
!
            DO is = 1, nsolid

              ! ... Assemble computational stencil
              CALL rnb(u,us(:,is),ijk)      
              CALL rnb(v,vs(:,is),ijk)
              CALL rnb(w,ws(:,is),ijk)
              CALL nb ( dens, rlk(:,is), ijk ) 
!
             IF (ctu > 0) THEN
                CALL ctu1_rnb(u,us(:,is),ijk)
                CALL ctu1_rnb(v,vs(:,is),ijk)
                CALL ctu1_rnb(w,ws(:,is),ijk)
                CALL ctu1_nb(dens,rlk(:,is),ijk)
                IF (ctu > 1) THEN
                  CALL ctu2_rnb(u,us(:,is),ijk)
                  CALL ctu2_rnb(v,vs(:,is),ijk)
                  CALL ctu2_rnb(w,ws(:,is),ijk)
                  CALL ctu2_nb(dens,rlk(:,is),ijk)
                  IF (ctu > 2) THEN
                    CALL ctu3_rnb(u,us(:,is),ijk)
                    CALL ctu3_rnb(v,vs(:,is),ijk)
                    CALL ctu3_rnb(w,ws(:,is),ijk)
                    CALL ctu3_nb(dens,rlk(:,is),ijk)
                  END IF
                END IF
              END IF
!
              ! ... Mask non-physical velocities from Immersed Boundaries
              !
              IF (immb >= 1) THEN
                IF (flag(ijk)==bl_cell .OR. flag(ijk)==immb_cell) THEN
                  u = maskstencil(u,mask_x)
                  v = maskstencil(v,mask_y)
                  w = maskstencil(w,mask_z)
                END IF
              END IF
!
              ! ... Interpolate density on the staggered grid
              CALL interpolate_x(dens, dens_stagx, i)
              CALL interpolate_y(dens, dens_stagy, j)
              CALL interpolate_z(dens, dens_stagz, k)

              CALL flu(usfe(ijk,is), usfn(ijk,is), usft(ijk,is),           &
                      usfe(imjk,is), usfn(ijmk,is), usft(ijkm,is),        &
                      dens_stagx, u, v, w, i)

              CALL flv(vsfe(ijk,is), vsfn(ijk,is), vsft(ijk,is),           &
                      vsfe(imjk,is), vsfn(ijmk,is), vsft(ijkm,is),        &
                      dens_stagy, u, v, w, j)

              CALL flw(wsfe(ijk,is), wsfn(ijk,is), wsft(ijk,is),           &
                      wsfe(imjk,is), wsfn(ijmk,is), wsft(ijkm,is),        &
                      dens_stagz, u, v, w, k)
!
! ... Second order MUSCL correction
!
              IF (muscl > 0 .AND. flag(ijk)==fluid) THEN
                IF ( i /= nx-1 ) &
                  CALL muscl_flu(usfe(ijk,is), usfn(ijk,is), usft(ijk,is), &
                           dens_stagx, u, v, w, i, j, k)
                IF ( j /= ny-1 ) &
                  CALL muscl_flv(vsfe(ijk,is), vsfn(ijk,is), vsft(ijk,is), &
                           dens_stagy, u, v, w, i, j, k)
                IF ( k /= nz-1 ) &
                  CALL muscl_flw(wsfe(ijk,is), wsfn(ijk,is), wsft(ijk,is), &
                           dens_stagz, u, v, w, i, j, k)
              END IF
!
              IF (ctu > 0) THEN
!
! ... First order Corner Transport Upwind correction
!
                CALL ctu1_flu(usfe(ijk,is), usfn(ijk,is), usft(ijk,is), &
                          dens_stagx, u, v, w, i, j, k)
                CALL ctu1_flv(vsfe(ijk,is), vsfn(ijk,is), vsft(ijk,is), &
                          dens_stagy, u, v, w, i, j, k)
                CALL ctu1_flw(wsfe(ijk,is), wsfn(ijk,is), wsft(ijk,is), &
                          dens_stagz, u, v, w, i, j, k)
!
! ... Second order Corner Transport Upwind correction
!
                IF (ctu > 1 .AND. muscl == 0 .AND. flag(ijk) == fluid) THEN
                  IF (i /= nx-1) CALL ctu2_flu(usfe(ijk,is), usfn(ijk,is), usft(ijk,is), &
                                           dens_stagx, u, v, w, i, j, k)
                  IF (j /= ny-1) CALL ctu2_flv(vsfe(ijk,is), vsfn(ijk,is), vsft(ijk,is), &
                                           dens_stagy, u, v, w, i, j, k)
                  IF (k /= nz-1) CALL ctu2_flw(wsfe(ijk,is), wsfn(ijk,is), wsft(ijk,is), &
                                           dens_stagz, u, v, w, i, j, k)
!
! ... Third order Corner Transport Upwind correction
!
                  IF (ctu > 2) THEN
                    CALL ctu3_flu(usfe(ijk,is), usfn(ijk,is), usft(ijk,is), &
                              dens_stagx, u, v, w, i, j, k)
                    CALL ctu3_flv(vsfe(ijk,is), vsfn(ijk,is), vsft(ijk,is), &
                              dens_stagy, u, v, w, i, j, k)
                    CALL ctu3_flw(wsfe(ijk,is), wsfn(ijk,is), wsft(ijk,is), &
                              dens_stagz, u, v, w, i, j, k)
                  END IF
                END IF
              END IF

            END DO

          END IF         
        END IF         

      END DO 
!
      CALL data_exchange(ugfe)
      CALL data_exchange(ugft)
      CALL data_exchange(wgfe)
      CALL data_exchange(wgft)

      CALL data_exchange(usfe)
      CALL data_exchange(usft)
      CALL data_exchange(wsfe)
      CALL data_exchange(wsft)

      IF (job_type == JOB_TYPE_3D) THEN
        CALL data_exchange(ugfn)
        CALL data_exchange(vgfe)
        CALL data_exchange(vgfn)
        CALL data_exchange(vgft)
        CALL data_exchange(wgfn)

        CALL data_exchange(usfn)
        CALL data_exchange(vsfe)
        CALL data_exchange(vsfn)
        CALL data_exchange(vsft)
        CALL data_exchange(wsfn)
      END IF
!
      RETURN
      END SUBROUTINE compute_all_fluxes_tilde
!----------------------------------------------------------------------
      SUBROUTINE test_fluxes_tilde
      USE domain_mapping, ONLY: ncint, meshinds
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm
      USE gas_solid_velocity, ONLY: ug
      USE io_files, ONLY: tempunit

      IMPLICIT NONE
      INTEGER :: ijk, imesh, i, j, k

      OPEN(UNIT=tempunit,FILE='pdac.fl',STATUS='UNKNOWN')
      WRITE(tempunit,*) 'Test momentum fluxes ...'
      WRITE(tempunit,*)
      DO ijk = 1, ncint
        CALL subscr(ijk)
        CALL meshinds(ijk,imesh,i,j,k)
          
        IF (i == 2 .AND. k == 31) THEN
          WRITE(tempunit,101) ugfe(ijk), ugft(ijk), ugfe(imjk), ugft(ijkm)
          WRITE(tempunit,101) wgfe(ijk), wgft(ijk), wgfe(imjk), wgft(ijkm)
          WRITE(tempunit,101)
        END IF
        IF (i == 1 .AND. k == 32) THEN
          WRITE(tempunit,101) ugfe(ijk), ugft(ijk), ugfe(imjk), ugft(ijkm)
          WRITE(tempunit,101) wgfe(ijk), wgft(ijk), wgfe(imjk), wgft(ijkm)
          WRITE(tempunit,101)
        END IF
        IF (i == 2 .AND. k == 32) THEN
          WRITE(tempunit,101) ugfe(ijk), ugft(ijk), ugfe(imjk), ugft(ijkm)
          WRITE(tempunit,101) wgfe(ijk), wgft(ijk), wgfe(imjk), wgft(ijkm)
          WRITE(tempunit,101)
        END IF

      END DO

 100  FORMAT(4(I5),3(G20.10E2))
 101  FORMAT(3(G20.10E2))
      CLOSE(tempunit)

      RETURN
      END SUBROUTINE test_fluxes_tilde
!----------------------------------------------------------------------
      END MODULE tilde_momentum
!----------------------------------------------------------------------
