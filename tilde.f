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
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_momentum
      USE dimensions
      USE domain_decomposition, ONLY: ncint
      USE control_flags, ONLY: job_type
      IMPLICIT NONE
!
      ALLOCATE( rugn(ncint), rwgn(ncint))
      ALLOCATE( rusn(ncint,nsolid), rwsn(ncint,nsolid))
      rugn = 0.D0; rwgn = 0.D0
      rusn = 0.D0; rwsn = 0.D0
      IF (job_type == '3D') THEN
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
      USE domain_decomposition, ONLY: ncint, ncdom
      USE control_flags, ONLY: job_type
      USE tilde_energy, ONLY: rhg, rhs
      IMPLICIT NONE
!
! ... Allocate and initialize "tilde" terms (gas)
!
      ALLOCATE(rug(ncdom), rwg(ncdom))
      rug = 0.0D0
      rwg = 0.0D0
      IF (job_type == '3D') THEN
        ALLOCATE(rvg(ncdom))
        rvg = 0.0D0
      END IF
      ALLOCATE(rhg(ncint))
      rhg = 0.D0
!
! ... Allocate and initialize "tilde" terms (particles).
!
      ALLOCATE(rus(ncdom,nsolid), rws(ncdom,nsolid))
      rus = 0.0D0
      rws = 0.0D0
      IF (job_type == '3D') THEN
        ALLOCATE(rvs(ncdom,nsolid))
        rvs = 0.0D0
      END IF
      ALLOCATE(rhs(ncint,nsolid))
      rhs = 0.D0
!
! ... Allocate and initialize interphase terms.
!
      ALLOCATE(appu(ncdom, ((nsolid+1)**2+(nsolid+1))/2),   &
               appw(ncdom, ((nsolid+1)**2+(nsolid+1))/2))

      appu = 0.0D0
      appw = 0.0D0
      IF (job_type == '3D') THEN
        ALLOCATE(appv(ncdom, ((nsolid+1)**2+(nsolid+1))/2) )
        appv = 0.0D0
      END IF
!
      RETURN
      END SUBROUTINE allocate_fluxes
!-----------------------------------------------------
      SUBROUTINE deallocate_fluxes

        USE control_flags, ONLY: job_type
        USE tilde_energy, ONLY: rhg, rhs
        
        DEALLOCATE(rug, rwg)
        DEALLOCATE(rus, rws)
        DEALLOCATE(appu, appw)
        IF (job_type == '3D') THEN
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
! ... Compute explicitly and store all fields and physical parameters
! ... at time ndt for explicit time-advancement
! ... (2D-3D-Compliant)
!
      USE control_flags, ONLY: job_type
      USE dimensions
      USE domain_decomposition, ONLY: ncint, myijk, ncdom, data_exchange
      USE domain_decomposition, ONLY: meshinds
      USE eos_gas, ONLY: rgpgc, rgpgcn, xgc
      USE gas_solid_density, ONLY: rgp, rgpn, rlk, rlkn
      USE gas_solid_temperature, ONLY: sieg, siegn, sies, siesn, tg
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_viscosity, ONLY: viscon, mug, kapg
      USE gas_solid_viscosity, ONLY: gas_viscosity, part_viscosity
      USE grid, ONLY: dz, dy, dx, fl_l
      USE grid, ONLY: indx, indy, indz
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE pressure_epsilon, ONLY: p, pn
      USE set_indexes
!
      IMPLICIT NONE
      SAVE
!
      REAL*8 :: dxp, dyp, dzp, indxp, indyp, indzp
      REAL*8 :: rgp_e, rgp_n, rgp_t, rlk_e, rlk_n, rlk_t
      INTEGER :: i, j, k, ijk, imesh, is, ig
!
      CALL data_exchange(rgp)
      CALL data_exchange(rlk)
!
      DO ijk = 1, ncint
!       IF(fl_l(ijk) == 1) THEN
          CALL meshinds(ijk,imesh,i,j,k)
          CALL subscr_fieldn(ijk)
!
          IF (job_type == '2D') THEN

            dxp=dx(i)+dx(i+1)
            dzp=dz(k)+dz(k+1)
            indxp=1.D0/dxp
            indzp=1.D0/dzp

            rgp_e = (dx(i+1)*rgp(ijk)+dx(i)*rgp(ijke))*indxp
            rgp_t = (dz(k+1)*rgp(ijk)+dz(k)*rgp(ijkt))*indzp
            rugn(ijk)  = rgp_e * ug(ijk)
            rwgn(ijk)  = rgp_t * wg(ijk)
           
            DO is = 1, nsolid
               rlk_e = (rlk(ijk,is)*dx(i+1)+rlk(ijke,is)*dx(i))*indxp
               rlk_t = (rlk(ijk,is)*dz(k+1)+rlk(ijkt,is)*dz(k))*indzp
               rusn(ijk,is)  = rlk_e * us(ijk,is)
               rwsn(ijk,is)  = rlk_t * ws(ijk,is)
            END DO
 
          ELSE IF (job_type == '3D') THEN

            dxp=dx(i)+dx(i+1)
            dyp=dy(j)+dy(j+1)
            dzp=dz(k)+dz(k+1)
            indxp=1.D0/dxp
            indyp=1.D0/dyp
            indzp=1.D0/dzp

            rgp_e = (dx(i+1)*rgp(ijk)+dx(i)*rgp(ijke))*indxp
            rgp_n = (dy(j+1)*rgp(ijk)+dy(j)*rgp(ijkn))*indyp
            rgp_t = (dz(k+1)*rgp(ijk)+dz(k)*rgp(ijkt))*indzp
            rugn(ijk)  = rgp_e * ug(ijk)
            rvgn(ijk)  = rgp_n * vg(ijk)
            rwgn(ijk)  = rgp_t * wg(ijk)
!     
            DO is = 1, nsolid
              rlk_e = (rlk(ijk,is)*dx(i+1)+rlk(ijke,is)*dx(i))*indxp
              rlk_n = (rlk(ijk,is)*dy(j+1)+rlk(ijkn,is)*dy(j))*indyp
              rlk_t = (rlk(ijk,is)*dz(k+1)+rlk(ijkt,is)*dz(k))*indzp
              rusn(ijk,is)  = rlk_e * us(ijk,is)
              rvsn(ijk,is)  = rlk_n * vs(ijk,is)
              rwsn(ijk,is)  = rlk_t * ws(ijk,is)
            END DO
         
          END IF 
!
          pn(ijk)    = p(ijk)
          rgpn(ijk)  = rgp(ijk)
          siegn(ijk) = sieg(ijk)
 
          DO ig=1,ngas
            rgpgcn(ijk,ig) = rgpgc(ijk,ig)
          END DO
!
          DO is = 1, nsolid
            rlkn(ijk,is)  = rlk(ijk,is)
            siesn(ijk,is) = sies(ijk,is)
          END DO

!        END IF

      END DO
!
! ... Compute the temperature-dependent gas viscosity and th. conductivity
!
      DO ijk = 1, ncint
         CALL viscon(mug(ijk), kapg(ijk), xgc(:,ijk), tg(ijk))
      END DO
!
      RETURN
      END SUBROUTINE fieldn
!----------------------------------------------------------------------
      SUBROUTINE tilde
!
      USE atmosphere, ONLY: gravz
      USE dimensions
      USE domain_decomposition, ONLY: meshinds
      USE domain_decomposition, ONLY: ncint, myijk, ncdom, data_exchange
      USE control_flags, ONLY: job_type
      USE gas_solid_density, ONLY: rog, rgp, rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE grid, ONLY: dx, dy, dz, fl_l
      USE grid, ONLY: indx, indy, indz, inx, inxb
      USE momentum_transfer, ONLY: kdrags, inter
      USE pressure_epsilon, ONLY: ep, p
      USE time_parameters, ONLY: dt, time
      USE turbulence_model, ONLY: iss, iturb
      USE gas_solid_viscosity, ONLY: viscg, viscs
      USE gas_solid_viscosity, ONLY: gas_viscosity, part_viscosity
      USE gas_solid_viscosity, ONLY: mug
      USE gas_solid_viscosity, ONLY: gvisx, gvisy, gvisz, pvisx, pvisy, pvisz
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: subscr, subscr_tilde, imjk, ijmk, ijkm, ijkt
!
      IMPLICIT NONE
      SAVE
!
      INTEGER :: i, j, k, is, imesh
      INTEGER :: ijk
      REAL*8 :: dxp, dyp, dzp, indxp, indyp, indzp
      REAL*8 :: ugfw, ugfs, ugfb, vgfw, vgfs, vgfb, wgfw, wgfs, wgfb
      REAL*8 :: ugfx, ugfy, ugfz, vgfx, vgfy, vgfz, wgfx, wgfy, wgfz
      REAL*8 :: usfw, usfs, usfb, vsfw, vsfs, vsfb, wsfw, wsfs, wsfb
      REAL*8 :: usfx, usfy, usfz, vsfx, vsfy, vsfz, wsfx, wsfy, wsfz
      REAL*8 :: dugs, dvgs, dwgs
      REAL*8 :: rug_tmp, rvg_tmp, rwg_tmp
      REAL*8 :: rus_tmp, rvs_tmp, rws_tmp
      REAL*8, ALLOCATABLE :: nul(:)
!
! ... Allocate and initialize gas and particle viscous stresses
!
      ALLOCATE(gvisx(ncint), gvisz(ncint))
      ALLOCATE(pvisx(ncint,nsolid), pvisz(ncint,nsolid))
      gvisx = 0.D0; gvisz = 0.D0
      pvisx = 0.D0; pvisz = 0.D0
      IF (job_type == '3D') THEN
        ALLOCATE(gvisy(ncint))
        ALLOCATE(pvisy(ncint,nsolid))
        gvisy = 0.D0
        pvisy = 0.D0
      END IF
!
      IF (iturb == 0) THEN
        CALL data_exchange(ug)
        CALL data_exchange(wg)
        IF (job_type == '3D') CALL data_exchange(vg)
      END IF
      IF (iss == 0) THEN
        CALL data_exchange(us)
        CALL data_exchange(ws)
        IF (job_type == '3D') CALL data_exchange(vs)
      END IF
      CALL data_exchange(ep)
!
! ... Calculate gas viscous stress tensor
!
      IF (gas_viscosity) CALL viscg        
!
! ... Calculate particles viscous stress tensor
!
      IF (part_viscosity) CALL viscs
!
! ... Allocate and initialize gas convective fluxes
!
      ALLOCATE(ugfe(ncdom), ugft(ncdom))
      ALLOCATE(wgfe(ncdom), wgft(ncdom))
      ugfe = 0.0D0; ugft = 0.0D0
      ugfw = 0.0D0; ugfb = 0.0D0
      wgfe = 0.0D0; wgft = 0.0D0
      wgfw = 0.0D0; wgfb = 0.0D0

      IF (job_type == '3D') THEN
        ALLOCATE(ugfn(ncdom))
        ALLOCATE(vgfe(ncdom), vgfn(ncdom), vgft(ncdom))
        ALLOCATE(wgfn(ncdom))
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
      ALLOCATE(usfe(ncdom,nsolid), usft(ncdom,nsolid))
      ALLOCATE(wsfe(ncdom,nsolid), wsft(ncdom,nsolid))

      usfe = 0.0D0; usft = 0.0D0
      usfw = 0.0D0; usfb = 0.0D0
      wsfe = 0.0D0; wsft = 0.0D0
      wsfw = 0.0D0; wsfb = 0.0D0

      IF (job_type == '3D') THEN
        ALLOCATE(usfn(ncdom,nsolid))
        ALLOCATE(vsfe(ncdom,nsolid), vsfn(ncdom,nsolid), vsft(ncdom,nsolid))
        ALLOCATE(wsfn(ncdom,nsolid))
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
      ALLOCATE(kpgv(nsolid))
      kpgv = 0.0D0
!
! ... (a temporary array used in 2D) ...
!
      ALLOCATE(nul( ((nsolid+1)**2+(nsolid+1))/2) )
      nul = 0.D0
!
! ... Compute all convective East, North, and Top fluxes 
! ... in the physical domain and ghost cells

      IF ( ( nsolid == 2 ) .AND. ( job_type == '3D' ) ) THEN   
        CALL compute_all_fluxes_3d_3phase
      ELSE
        CALL compute_all_fluxes
      END IF
!
! ... Fluxes on West, South and Bottom sides keep values 
! ... of East, North and Top fluxes from neighbouring cells.
!
      DO ijk = 1, ncint
        IF(fl_l(ijk) == 1) THEN
          CALL meshinds(ijk,imesh,i,j,k)
          CALL subscr_tilde(ijk)
!
          ugfw = ugfe(imjk)
          ugfb = ugft(ijkm)
          wgfw = wgfe(imjk)
          wgfb = wgft(ijkm)
!
          ugfx = ugfe(ijk) - ugfw
          ugfz = ugft(ijk) - ugfb
          wgfx = wgfe(ijk) - wgfw
          wgfz = wgft(ijk) - wgfb
!
          IF (job_type == '2D') THEN

            ugfy = 0.D0
            vgfx = 0.D0
            vgfy = 0.D0
            vgfz = 0.D0
            wgfy = 0.D0

          ELSE IF (job_type == '3D') THEN

            ugfs = ugfn(ijmk)
            vgfw = vgfe(imjk)
            vgfs = vgfn(ijmk)
            vgfb = vgft(ijkm)
            wgfs = wgfn(ijmk)
!
            ugfy = ugfn(ijk) - ugfs
            vgfx = vgfe(ijk) - vgfw
            vgfy = vgfn(ijk) - vgfs
            vgfz = vgft(ijk) - vgfb
            wgfy = wgfn(ijk) - wgfs
            
          END IF
!
! ... compute explicit (tilde) terms in the momentum equation (gas)
! 
          dxp=dx(i)+dx(i+1)
          dyp=dy(j)+dy(j+1)
          dzp=dz(k)+dz(k+1)
          indxp=1.D0/dxp
          indzp=1.D0/dzp
          indyp=1.D0/dyp
!         
          rug_tmp = rugn(ijk) + dt * gvisx(ijk)                     
          rug_tmp = rug_tmp - dt * indxp * 2.D0 * ugfx * inxb(i)        
          rug_tmp = rug_tmp - dt * indy(j) * ugfy                  
          rug_tmp = rug_tmp - dt * indz(k) * ugfz   
          rug (ijk) = rug_tmp
!
          rwg_tmp = rwgn(ijk) + dt * gvisz(ijk)                     
          rwg_tmp = rwg_tmp + dt * (dz(k+1)*rgp(ijk)+dz(k)*rgp(ijkt))*indzp * gravz  
          rwg_tmp = rwg_tmp - dt * indx(i) * wgfx * inx(i)                           
          rwg_tmp = rwg_tmp - dt * indy(j) * wgfy                                    
          rwg_tmp = rwg_tmp - dt * indzp * 2.D0 * wgfz
          rwg(ijk) = rwg_tmp
!
          IF (job_type == '3D') THEN

            rvg_tmp = rvgn(ijk) + dt * gvisy(ijk)                     
            rvg_tmp = rvg_tmp - dt * indx(i) * vgfx               
            rvg_tmp = rvg_tmp - dt * indyp * 2.D0 * vgfy   
            rvg_tmp = rvg_tmp - dt * indz(k) * vgfz    
            rvg(ijk) = rvg_tmp

          END IF

! ... same procedure carried out for particulate phases


          IF (job_type == '2D') THEN

             DO is = 1, nsolid
!
! ... West, South and Bottom fluxes (particles)
!
                usfw = usfe(imjk,is)
                usfx = usfe(ijk,is) - usfw
                usfb = usft(ijkm,is)
                usfz = usft(ijk,is) - usfb
                wsfw = wsfe(imjk,is)
                wsfx = wsfe(ijk,is) - wsfw
                wsfb = wsft(ijkm,is)
                wsfz = wsft(ijk,is) - wsfb
                usfy = 0.D0
                vsfx = 0.D0
                vsfy = 0.D0
                vsfz = 0.D0
                wsfy = 0.D0
!
! ... compute explicit (tilde) terms in the momentum equation (particles)
! 
                rus_tmp = rusn(ijk,is) + dt*pvisx(ijk,is)               
                rus_tmp = rus_tmp - dt*indxp*2.D0* usfx * inxb(i)   
                rus_tmp = rus_tmp - dt*indy(j)* usfy  
                rus_tmp = rus_tmp - dt*indz(k)* usfz 
                rus(ijk,is) = rus_tmp                   
!
                rws_tmp = rwsn(ijk,is) + dt*pvisz(ijk,is)              
                rws_tmp = rws_tmp + dt*(rlk(ijk,is)*dz(k+1)+rlk(ijkt,is)*dz(k))*indzp*gravz 
                rws_tmp = rws_tmp - dt*indx(i)* wsfx * inx(i)                 
                rws_tmp = rws_tmp - dt*indy(j)* wsfy                      
                rws_tmp = rws_tmp - dt*indzp*2.D0* wsfz  
                rws(ijk,is) = rws_tmp
!
! ... Compute the gas-particle drag coefficients
!
                dugs = ( (ug(ijk)-us(ijk,is)) + (ug(imjk)-us(imjk,is)) )*0.5D0
                dwgs = ( (wg(ijk)-ws(ijk,is)) + (wg(ijkm)-ws(ijkm,is)) )*0.5D0
                dvgs = 0.D0
                CALL kdrags(kpgv(is), dugs, dvgs, dwgs, ep(ijk),         &
                        rgp(ijk), rlk(ijk,is), mug(ijk), is)                  
            END DO 

          ELSE IF (job_type == '3D') THEN

            DO is = 1, nsolid
!
! ... West, South and Bottom fluxes (particles)
!
               usfw = usfe(imjk,is)
               usfx = usfe(ijk,is) - usfw
               usfs = usfn(ijmk,is)
               usfy = usfn(ijk,is) - usfs
               usfb = usft(ijkm,is)
               usfz = usft(ijk,is) - usfb
               vsfw = vsfe(imjk,is)
               vsfx = vsfe(ijk,is) - vsfw
               vsfs = vsfn(ijmk,is)
               vsfy = vsfn(ijk,is) - vsfs
               vsfb = vsft(ijkm,is)
               vsfz = vsft(ijk,is) - vsfb
               wsfw = wsfe(imjk,is)
               wsfx = wsfe(ijk,is) - wsfw
               wsfb = wsft(ijkm,is)
               wsfz = wsft(ijk,is) - wsfb 
               wsfs = wsfn(ijmk,is)             
               wsfy = wsfn(ijk,is) - wsfs
!
! ... compute explicit (tilde) terms in the momentum equation (particles)
! 
               rus_tmp = rusn(ijk,is) + dt*pvisx(ijk,is)               
               rus_tmp = rus_tmp - dt*indxp*2.D0* usfx * inxb(i)   
               rus_tmp = rus_tmp - dt*indy(j)* usfy  
               rus_tmp = rus_tmp - dt*indz(k)* usfz 
               rus(ijk,is) = rus_tmp                   
!
               rvs_tmp = rvsn(ijk,is) + dt*pvisy(ijk,is)              
               rvs_tmp = rvs_tmp - dt*indx(i)* vsfx       
               rvs_tmp = rvs_tmp - dt*indyp*2.D0* vsfy              
               rvs_tmp = rvs_tmp - dt*indz(k)* vsfz    
               rvs(ijk,is) = rvs_tmp
!
               rws_tmp = rwsn(ijk,is) + dt*pvisz(ijk,is)              
               rws_tmp = rws_tmp + dt*(rlk(ijk,is)*dz(k+1)+rlk(ijkt,is)*dz(k))*indzp*gravz
               rws_tmp = rws_tmp - dt*indx(i)* wsfx * inx(i)                 
               rws_tmp = rws_tmp - dt*indy(j)* wsfy                      
               rws_tmp = rws_tmp - dt*indzp*2.D0* wsfz  
               rws(ijk,is) = rws_tmp
!
! ... Compute the gas-particle drag coefficients
!
               dugs = ( (ug(ijk)-us(ijk,is)) + (ug(imjk)-us(imjk,is)) )*0.5D0
               dvgs = ( (vg(ijk)-vs(ijk,is)) + (vg(ijmk)-vs(ijmk,is)) )*0.5D0
               dwgs = ( (wg(ijk)-ws(ijk,is)) + (wg(ijkm)-ws(ijkm,is)) )*0.5D0  
           
   
               CALL kdrags(kpgv(is), dugs, dvgs, dwgs, ep(ijk),         &
                        rgp(ijk), rlk(ijk,is), mug(ijk), is)                  
            END DO

          END IF

 
!
! ... Compute the particle-particle coefficients and the interphase matrix
!
          IF (job_type == '2D') THEN
            CALL inter(appu(ijk,:), nul(:), appw(ijk,:), kpgv(:),    &
     &                 us, us, ws, rlk, ijk)
          ELSE IF (job_type == '3D') THEN
            CALL inter(appu(ijk,:), appv(ijk,:), appw(ijk,:), kpgv(:),    &
     &                 us, vs, ws, rlk, ijk)
          END IF
!
        END IF
      END DO
!
!      CALL test_fluxes
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
      DEALLOCATE(nul)
!
      CALL data_exchange(appu)
      CALL data_exchange(appw)
      CALL data_exchange(rug)
      CALL data_exchange(rwg)
      CALL data_exchange(rus)
      CALL data_exchange(rws)
!
      IF (job_type == '3D') THEN

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
      SUBROUTINE compute_all_fluxes
!
      USE atmosphere, ONLY: gravz
      USE dimensions
      USE domain_decomposition, ONLY: meshinds
      USE domain_decomposition, ONLY: ncint, myijk, data_exchange
      USE convective_fluxes_u, ONLY: flu
      USE convective_fluxes_v, ONLY: flv
      USE convective_fluxes_w, ONLY: flw
      USE control_flags, ONLY: job_type
      USE gas_solid_density, ONLY: rog, rgp, rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE grid, ONLY: dx, dy, dz, fl_l
      USE grid, ONLY: indx, indy, indz, inx, inxb
      USE pressure_epsilon, ONLY: ep, p
      USE time_parameters, ONLY: dt, time
      USE turbulence_model, ONLY: iss, iturb
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm, ijkt, imjpk, imjkp, ipjmk, ijmkp, ipjkm, ijpkm 
      USE set_indexes, ONLY: stencil, nb, rnb, rnb_13
!
      IMPLICIT NONE
      SAVE
!
      INTEGER :: i, j, k, is, imesh
      INTEGER :: ijk
      REAL*8 :: dxp, dyp, dzp, indxp, indyp, indzp
      TYPE(stencil) :: u, v, w, dens
!
! ... Compute fluxes on East, North and Top sides of a cell
! ... in the whole computational domain.
!
      DO ijk = 1, ncint
        IF(fl_l(ijk) == 1) THEN
          CALL subscr(ijk)
          CALL meshinds(ijk, imesh,i,j,k)
!
          IF (job_type == '2D') THEN

            CALL nb(dens,rgp,ijk)
            CALL rnb(u,ug,ijk)
            CALL rnb(w,wg,ijk)

            CALL flu(ugfe(ijk), ugft(ijk),                  &
                     ugfe(imjk), ugft(ijkm), dens, u, w, ijk)
            CALL flw(wgfe(ijk), wgft(ijk),                   &
                      wgfe(imjk), wgft(ijkm), dens, u, w, ijk)

            DO is = 1, nsolid
              CALL nb(dens,rlk(:,is),ijk)
              CALL rnb(u,us(:,is),ijk)
              CALL rnb(w,ws(:,is),ijk)

              CALL flu(usfe(ijk,is), usft(ijk,is),                     &
                       usfe(imjk,is), usft(ijkm,is), dens, u, w, ijk)
              CALL flw(wsfe(ijk,is), wsft(ijk,is),                     &
                       wsfe(imjk,is), wsft(ijkm,is), dens, u, w, ijk)
            END DO
            
          ELSE IF (job_type == '3D') THEN

            CALL nb ( dens, rgp, ijk )
            !CALL rnb(u,ug,ijk)
            CALL rnb_13 ( u, ug, ijk )
            u%wn = ug ( imjpk ) 
            u%wt = ug ( imjkp )
            !CALL rnb(w,wg,ijk) 
            CALL rnb_13( v, vg, ijk )
            v%es = vg ( ipjmk )
            v%st = vg ( ijmkp )
            !CALL rnb(v,vg,ijk)
            CALL rnb_13( w, wg, ijk )
            w%eb = wg ( ipjkm )
            w%nb = wg ( ijpkm )

            CALL flu(ugfe(ijk), ugfn(ijk), ugft(ijk),                     &
                    ugfe(imjk), ugfn(ijmk), ugft(ijkm), dens, u, v, w, ijk)
            CALL flv(vgfe(ijk), vgfn(ijk), vgft(ijk),                        &
                    vgfe(imjk), vgfn(ijmk), vgft(ijkm), dens, u, v, w, ijk)
            CALL flw(wgfe(ijk), wgfn(ijk), wgft(ijk),                        &
                    wgfe(imjk), wgfn(ijmk), wgft(ijkm), dens, u, v, w, ijk)
!
            DO is = 1, nsolid

              CALL nb ( dens, rlk(:,is), ijk ) 
              !CALL rnb(u,us(:,is),ijk)      
              CALL rnb_13 ( u, us(:,is), ijk )
              u%wn = us ( imjpk, is ) 
              u%wt = us ( imjkp, is) 
              !CALL rnb(v,vs(:,is),ijk)
              CALL rnb_13( v, vs(:,is), ijk )
              v%es = vs ( ipjmk, is )
              v%st = vs ( ijmkp, is )
              !CALL rnb(w,ws(:,is),ijk)
              CALL rnb_13( w, ws(:,is) , ijk )
              w%eb = ws ( ipjkm, is )
              w%nb = ws ( ijpkm, is )
!
              CALL flu(usfe(ijk,is), usfn(ijk,is), usft(ijk,is),           &
                      usfe(imjk,is), usfn(ijmk,is), usft(ijkm,is),        &
                      dens, u, v, w, ijk)
              CALL flv(vsfe(ijk,is), vsfn(ijk,is), vsft(ijk,is),           &
                      vsfe(imjk,is), vsfn(ijmk,is), vsft(ijkm,is),        &
                      dens, u, v, w, ijk)
              CALL flw(wsfe(ijk,is), wsfn(ijk,is), wsft(ijk,is),           &
                      wsfe(imjk,is), wsfn(ijmk,is), wsft(ijkm,is),        &
                      dens, u, v, w, ijk)
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

      IF (job_type == '3D') THEN
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
      END SUBROUTINE compute_all_fluxes

!----------------------------------------------------------------------
      SUBROUTINE compute_all_fluxes_3d_3phase
!
      USE atmosphere, ONLY: gravz
      USE dimensions
      USE domain_decomposition, ONLY: meshinds
      USE domain_decomposition, ONLY: ncint, myijk, data_exchange
      USE convective_fluxes_u, ONLY: flu, flu_1st
      USE convective_fluxes_v, ONLY: flv, flv_1st
      USE convective_fluxes_w, ONLY: flw, flw_1st
      USE flux_limiters, ONLY: muscl
      USE gas_solid_density, ONLY: rog, rgp, rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE grid, ONLY: dx, dy, dz, fl_l
      USE grid, ONLY: indx, indy, indz, inx, inxb
      USE pressure_epsilon, ONLY: ep, p
      USE time_parameters, ONLY: dt, time
      USE turbulence_model, ONLY: iss, iturb
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes
!
      IMPLICIT NONE
      SAVE
!
      INTEGER :: i, j, k, is, imesh
      INTEGER :: ijk
      REAL*8 :: dxp, dyp, dzp, indxp, indyp, indzp
      TYPE(stencil) :: u, v, w, dens
      TYPE(stencil) :: u1, v1, w1, dens1
      TYPE(stencil) :: u2, v2, w2, dens2
!
! ... Compute fluxes on East, North and Top sides of a cell
! ... in the whole computational domain.
!
      DO ijk = 1, ncint

        IF(fl_l(ijk) == 1) THEN

          CALL subscr(ijk)
          CALL meshinds(ijk, imesh,i,j,k)
!
          IF ( muscl == 0) THEN

            !CALL nb ( dens, rgp, ijk )
            dens%c = rgp ( ijk ) 
            dens%w = rgp ( ijkw ) 
            dens%s = rgp ( ijks ) 
            dens%b = rgp ( ijkb ) 
            dens%e = rgp ( ijke ) 
            dens%n = rgp ( ijkn ) 
            dens%t = rgp ( ijkt )  
            dens%ee = rgp ( ijkee ) 
            dens%en = rgp ( ijken ) 
            dens%et = rgp ( ijket ) 
            dens%es = rgp ( ijkes ) 
            dens%eb = rgp ( ijkeb ) 
            dens%nn = rgp ( ijknn ) 
            dens%nb = rgp ( ijknb ) 
            dens%nt = rgp ( ijknt ) 
            dens%wn = rgp ( ijkwn ) 
            dens%wt = rgp ( ijkwt ) 
            dens%tt = rgp ( ijktt ) 
            dens%st = rgp ( ijkst ) 

            !CALL rnb(u,ug,ijk)
            u%c = ug ( ijk ) 
            u%w = ug ( imjk ) 
            u%s = ug ( ijmk ) 
            u%b = ug ( ijkm ) 
            u%e = ug ( ipjk ) 
            u%n = ug ( ijpk ) 
            u%t = ug ( ijkp )    
            u%wn = ug ( imjpk ) 
            u%wt = ug ( imjkp )
          
            !CALL rnb(v,vg,ijk)
            v%c = vg ( ijk )
            v%w = vg ( imjk ) 
            v%s = vg ( ijmk ) 
            v%b = vg ( ijkm ) 
            v%e = vg ( ipjk ) 
            v%n = vg ( ijpk ) 
            v%t = vg ( ijkp ) 
            v%es = vg ( ipjmk )
            v%st = vg ( ijmkp )         
          
            !CALL rnb(w,wg,ijk) 
            w%c = wg ( ijk )
            w%w = wg ( imjk ) 
            w%s = wg ( ijmk ) 
            w%b = wg ( ijkm ) 
            w%e = wg ( ipjk ) 
            w%n = wg ( ijpk ) 
            w%t = wg ( ijkp ) 
            w%eb = wg ( ipjkm )
            w%nb = wg ( ijpkm )

            CALL flu_1st(ugfe(ijk), ugfn(ijk), ugft(ijk), ugfe(imjk),         &
                   ugfn(ijmk), ugft(ijkm), dens, u, v, w, i)
            CALL flv_1st(vgfe(ijk), vgfn(ijk), vgft(ijk), vgfe(imjk),         &
                    vgfn(ijmk), vgft(ijkm), dens, u, v, w, j)
            CALL flw_1st(wgfe(ijk), wgfn(ijk), wgft(ijk),                     &
                    wgfe(imjk), wgfn(ijmk), wgft(ijkm), dens, u, v, w, k)

          ELSE

            CALL nb ( dens, rgp, ijk )

            CALL first_rnb(u,ug,ijk)
            u%nn = ug ( ijppk ) 
            u%tt = ug ( ijkpp )
            u%ee = ug ( ippjk )
            u%wn = ug ( imjpk ) 
            u%wt = ug ( imjkp )

            CALL first_rnb(v,vg,ijk)
            v%nn = vg ( ijppk ) 
            v%tt = vg ( ijkpp )
            v%ee = vg ( ippjk )
            v%es = vg ( ipjmk )
            v%st = vg ( ijmkp ) 

            CALL first_rnb(w,wg,ijk) 
            w%nn = wg ( ijppk ) 
            w%tt = wg ( ijkpp )
            w%ee = wg ( ippjk )
            w%eb = wg ( ipjkm )
            w%nb = wg ( ijpkm )

            CALL flu(ugfe(ijk), ugfn(ijk), ugft(ijk), ugfe(imjk),             &
                    ugfn(ijmk), ugft(ijkm), dens, u, v, w, ijk )
            CALL flv(vgfe(ijk), vgfn(ijk), vgft(ijk), vgfe(imjk),             &
                    vgfn(ijmk), vgft(ijkm), dens, u, v, w, ijk)
            CALL flw(wgfe(ijk), wgfn(ijk), wgft(ijk),                        &
                    wgfe(imjk), wgfn(ijmk), wgft(ijkm), dens, u, v, w, ijk)
          END IF
!
          !DO is = 1, nsolid
!         Unrolled loop
         
          IF ( muscl == 0) THEN

            !CALL nb ( dens1, rlk(:,1), ijk ) 
            !CALL nb ( dens2, rlk(:,2), ijk )       
            dens1%c = rlk ( ijk,1 ) 
            dens1%w = rlk ( ijkw,1 ) 
            dens1%s = rlk ( ijks,1 )
            dens1%b = rlk ( ijkb,1 ) 
            dens1%e = rlk ( ijke,1 )
            dens1%n = rlk ( ijkn,1 ) 
            dens1%t = rlk ( ijkt,1 )
            dens1%ee = rlk ( ijkee,1 ) 
            dens1%en = rlk ( ijken,1 ) 
            dens1%et = rlk ( ijket,1 )
            dens1%es = rlk ( ijkes,1 ) 
            dens1%eb = rlk ( ijkeb,1 ) 
            dens1%nn = rlk ( ijknn,1 ) 
            dens1%nb = rlk ( ijknb,1 )
            dens1%nt = rlk ( ijknt,1 )
            dens1%wn = rlk ( ijkwn,1 ) 
            dens1%wt = rlk ( ijkwt,1 ) 
            dens1%tt = rlk ( ijktt,1 )         
            dens1%st = rlk ( ijkst,1 )  

            dens2%c = rlk ( ijk,2 ) 
            dens2%w = rlk ( ijkw,2 ) 
            dens2%s = rlk ( ijks,2 )
            dens2%b = rlk ( ijkb,2 ) 
            dens2%e = rlk ( ijke,2 ) 
            dens2%n = rlk ( ijkn,2 ) 
            dens2%t = rlk ( ijkt,2 )
            dens2%ee = rlk ( ijkee,2 ) 
            dens2%en = rlk ( ijken,2 ) 
            dens2%et = rlk ( ijket,2 ) 
            dens2%es = rlk ( ijkes,2 ) 
            dens2%eb = rlk ( ijkeb,2 )
            dens2%nn = rlk ( ijknn,2 )
            dens2%nb = rlk ( ijknb,2 ) 
            dens2%nt = rlk ( ijknt,2 )
            dens2%wn = rlk ( ijkwn,2 )
            dens2%wt = rlk ( ijkwt,2 ) 
            dens2%tt = rlk ( ijktt,2 )          
            dens2%st = rlk ( ijkst,2 )

            !CALL rnb_13 ( u1, us(:,1), ijk )
            !CALL rnb_13 ( u2, us(:,2), ijk ) 
            u1%c = us ( ijk,1 )   
            u1%w = us ( imjk,1 ) 
            u1%s = us ( ijmk,1 ) 
            u1%b = us ( ijkm,1 ) 
            u1%e = us ( ipjk,1 )
            u1%n = us ( ijpk,1 ) 
            u1%t = us ( ijkp,1 ) 
            u1%wn = us ( imjpk,1 ) 
            u1%wt = us ( imjkp,1 )

            u2%c = us ( ijk,2 ) 
            u2%w = us ( imjk,2 ) 
            u2%w = us ( imjk,2 ) 
            u2%s = us ( ijmk,2 ) 
            u2%b = us ( ijkm,2 ) 
            u2%e = us ( ipjk,2 ) 
            u2%n = us ( ijpk,2 ) 
            u2%t = us ( ijkp,2 ) 
            u2%wn = us ( imjpk,2 ) 
            u2%wt = us ( imjkp,2 )


            !CALL rnb_13( v1, vs(:,1), ijk )
            !CALL rnb_13( v2, vs(:,2), ijk )
            v1%c = vs ( ijk,1 )      
            v1%w = vs ( imjk,1 ) 
            v1%s = vs ( ijmk,1 ) 
            v1%b = vs ( ijkm,1 ) 
            v1%e = vs ( ipjk,1 )
            v1%n = vs ( ijpk,1 ) 
            v1%t = vs ( ijkp,1 ) 
            v1%es = vs ( ipjmk, 1 )
            v1%st = vs ( ijmkp, 1 )
      
            v2%c = vs ( ijk,2 )
            v2%w = vs ( imjk,2 ) 
            v2%s = vs ( ijmk,2 ) 
            v2%b = vs ( ijkm,2 ) 
            v2%e = vs ( ipjk,2 ) 
            v2%n = vs ( ijpk,2 ) 
            v2%t = vs ( ijkp,2 ) 
            v2%es = vs ( ipjmk, 2 )
            v2%st = vs ( ijmkp, 2 ) 
            
            !CALL rnb_13( w1, ws(:,1) , ijk )             
            !CALL rnb_13( w2, ws(:,2) , ijk )
            w1%c = ws ( ijk,1 )
            w1%w = ws ( imjk,1 ) 
            w1%s = ws ( ijmk,1 ) 
            w1%b = ws ( ijkm,1 ) 
            w1%e = ws ( ipjk,1 )
            w1%n = ws ( ijpk,1 ) 
            w1%t = ws ( ijkp,1 ) 
            w1%eb = ws ( ipjkm, 1 )
            w1%nb = ws ( ijpkm, 1 ) 
          
            w2%c = ws ( ijk,2 )
            w2%w = ws ( imjk,2 ) 
            w2%s = ws ( ijmk,2 ) 
            w2%b = ws ( ijkm,2 ) 
            w2%e = ws ( ipjk,2 ) 
            w2%n = ws ( ijpk,2 ) 
            w2%t = ws ( ijkp,2 )
            w2%eb = ws ( ipjkm, 2 )
            w2%nb = ws ( ijpkm, 2 ) 

            CALL flu_1st (usfe(ijk,1), usfn(ijk,1), usft(ijk,1),    &
                     usfe(imjk,1), usfn(ijmk,1), usft(ijkm,1),        &
                     dens1, u1, v1, w1, i)
            CALL flu_1st (usfe(ijk,2), usfn(ijk,2), usft(ijk,2),    &
                     usfe(imjk,2), usfn(ijmk,2), usft(ijkm,2),        &
                     dens2, u2, v2, w2, i)
            
            CALL flv_1st(vsfe(ijk,1), vsfn(ijk,1), vsft(ijk,1),     &
                     vsfe(imjk,1), vsfn(ijmk,1), vsft(ijkm,1),        &
                     dens1, u1, v1, w1, j)
            CALL flv_1st(vsfe(ijk,2), vsfn(ijk,2), vsft(ijk,2),     &
                     vsfe(imjk,2), vsfn(ijmk,2), vsft(ijkm,2),        &
                     dens2, u2, v2, w2, j)
            
            CALL flw_1st(wsfe(ijk,1), wsfn(ijk,1), wsft(ijk,1),     &
                     wsfe(imjk,1), wsfn(ijmk,1), wsft(ijkm,1),        &
                     dens1, u1, v1, w1, k)   
            CALL flw_1st(wsfe(ijk,2), wsfn(ijk,2), wsft(ijk,2),     &
                     wsfe(imjk,2), wsfn(ijmk,2), wsft(ijkm,2),        &
                     dens2, u2, v2, w2, k)
          ELSE

            CALL nb ( dens1, rlk(:,1), ijk ) 
            CALL nb ( dens2, rlk(:,2), ijk )       
         
            CALL first_rnb ( u1, us(:,1), ijk )
            CALL first_rnb ( u2, us(:,2), ijk ) 
            u1%nn = us ( ijppk,1 )
            u1%tt = us ( ijkpp,1 )
            u1%ee = us ( ippjk,1 )
            u1%wn = us ( imjpk,1 ) 
            u1%wt = us ( imjkp,1 )
           
            u2%nn = us ( ijppk,2 ) 
            u2%tt = us ( ijkpp,2 )
            u2%ee = us ( ippjk,2 )
            u2%wn = us ( imjpk,2 ) 
            u2%wt = us ( imjkp,2 )

            CALL first_rnb( v1, vs(:,1), ijk )
            CALL first_rnb( v2, vs(:,2), ijk )
            v1%nn = vs ( ijppk,1 )
            v1%tt = vs ( ijkpp,1 )
            v1%ee = vs ( ippjk,1 )
            v1%es = vs ( ipjmk,1 )
            v1%st = vs ( ijmkp,1 )

            v2%nn = vs ( ijppk,2 )
            v2%tt = vs ( ijkpp,2 )
            v2%ee = vs ( ippjk,2 )
            v2%es = vs ( ipjmk,2 )
            v2%st = vs ( ijmkp,2 )    

            CALL first_rnb( w1, ws(:,1) , ijk )             
            CALL first_rnb( w2, ws(:,2) , ijk )
            w1%nn = ws ( ijppk,1 )
            w1%tt = ws ( ijkpp,1 )
            w1%ee = ws ( ippjk,1 )         
            w1%eb = ws ( ipjkm,1 )
            w1%nb = ws ( ijpkm,1 ) 
          
            w2%nn = ws ( ijppk,2 )
            w2%tt = ws ( ijkpp,2 )
            w2%ee = ws ( ippjk,2 )
            w2%eb = ws ( ipjkm,2 )
            w2%nb = ws ( ijpkm,2 ) 

            CALL flu(usfe(ijk,1), usfn(ijk,1), usft(ijk,1),         &
                     usfe(imjk,1), usfn(ijmk,1), usft(ijkm,1),        &
                     dens1, u1, v1, w1, ijk)
            CALL flu(usfe(ijk,2), usfn(ijk,2), usft(ijk,2),         &
                     usfe(imjk,2), usfn(ijmk,2), usft(ijkm,2),        &
                    dens2, u2, v2, w2, ijk)
            
            CALL flv(vsfe(ijk,1), vsfn(ijk,1), vsft(ijk,1),         &
                     vsfe(imjk,1), vsfn(ijmk,1), vsft(ijkm,1),        &
                     dens1, u1, v1, w1, ijk)
            CALL flv(vsfe(ijk,2), vsfn(ijk,2), vsft(ijk,2),         &
                     vsfe(imjk,2), vsfn(ijmk,2), vsft(ijkm,2),        &
                     dens2, u2, v2, w2, ijk)
            
            CALL flw(wsfe(ijk,1), wsfn(ijk,1), wsft(ijk,1),         &
                     wsfe(imjk,1), wsfn(ijmk,1), wsft(ijkm,1),        &
                     dens1, u1, v1, w1, ijk)   
            CALL flw(wsfe(ijk,2), wsfn(ijk,2), wsft(ijk,2),         &
                     wsfe(imjk,2), wsfn(ijmk,2), wsft(ijkm,2),        &
                     dens2, u2, v2, w2, ijk)

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
    
!
      RETURN
      END SUBROUTINE compute_all_fluxes_3d_3phase
!----------------------------------------------------------------------

      SUBROUTINE test_fluxes
!
      USE control_flags, ONLY: job_type
      USE dimensions, ONLY: nsolid
      USE io_restart, ONLY: write_array
      USE kinds
      USE parallel, ONLY: nproc, mpime, root, group
      USE time_parameters, ONLY: time
!
      IMPLICIT NONE
!
      CHARACTER :: filnam*11
!
      INTEGER :: ig,is
      LOGICAL :: lform = .TRUE.
!
      filnam='output.test'

      IF( mpime == root ) THEN

        IF (lform) THEN
          OPEN(UNIT=12,FILE=filnam)
          WRITE(12,*) time
        ELSE 
          OPEN(UNIT=12,FORM='UNFORMATTED',FILE=filnam)
          WRITE(12) REAL(time,4)
        END IF
 
      END IF
!
      IF (job_type == '2D') THEN
        CALL write_array( 12, rug, sgl, lform )
        CALL write_array( 12, rwg, sgl, lform )
        CALL write_array( 12, ugfe, sgl, lform )
        CALL write_array( 12, ugft, sgl, lform )
        CALL write_array( 12, wgfe, sgl, lform )
        CALL write_array( 12, wgft, sgl, lform )
      ELSE IF (job_type == '3D') THEN
        CALL write_array( 12, rug, sgl, lform )
        CALL write_array( 12, rvg, sgl, lform )
        CALL write_array( 12, rwg, sgl, lform )
        CALL write_array( 12, ugfe, sgl, lform )
        CALL write_array( 12, ugfn, sgl, lform )
        CALL write_array( 12, ugft, sgl, lform )
        CALL write_array( 12, vgfe, sgl, lform )
        CALL write_array( 12, vgfn, sgl, lform )
        CALL write_array( 12, vgft, sgl, lform )
        CALL write_array( 12, wgfe, sgl, lform )
        CALL write_array( 12, wgfn, sgl, lform )
        CALL write_array( 12, wgft, sgl, lform )
      END IF
!
      DO is = 1, nsolid
        IF (job_type == '2D') THEN
          CALL write_array( 12, rus(:,is), sgl, lform )
          CALL write_array( 12, rws(:,is), sgl, lform )
          CALL write_array( 12, usfe(:,is), sgl, lform )
          CALL write_array( 12, usft(:,is), sgl, lform )
          CALL write_array( 12, wsfe(:,is), sgl, lform )
          CALL write_array( 12, wsft(:,is), sgl, lform )
        ELSE IF (job_type == '3D') THEN
          CALL write_array( 12, rus(:,is), sgl, lform )
          CALL write_array( 12, rvs(:,is), sgl, lform )
          CALL write_array( 12, rws(:,is), sgl, lform )
          CALL write_array( 12, usfe(:,is), sgl, lform )
          CALL write_array( 12, usfn(:,is), sgl, lform )
          CALL write_array( 12, usft(:,is), sgl, lform )
          CALL write_array( 12, vsfe(:,is), sgl, lform )
          CALL write_array( 12, vsfn(:,is), sgl, lform )
          CALL write_array( 12, vsft(:,is), sgl, lform )
          CALL write_array( 12, wsfe(:,is), sgl, lform )
          CALL write_array( 12, wsfn(:,is), sgl, lform )
          CALL write_array( 12, wsft(:,is), sgl, lform )
        END IF
      END DO

      IF( mpime == root ) THEN
        CLOSE (12)
      END IF
!
      RETURN
      END SUBROUTINE test_fluxes
!----------------------------------------------------------------------
      END MODULE tilde_momentum
!----------------------------------------------------------------------
