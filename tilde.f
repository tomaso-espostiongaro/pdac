!----------------------------------------------------------------------
      MODULE tilde_momentum
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: rug, rvg, rwg
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rus, rvs, rws
      REAL*8, DIMENSION(:),   ALLOCATABLE :: rugn, rvgn, rwgn
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rusn, rvsn, rwsn
!
      REAL*8, ALLOCATABLE :: ugfe(:), ugfn(:), ugft(:)
      REAL*8, ALLOCATABLE :: vgfe(:), vgfn(:), vgft(:)
      REAL*8, ALLOCATABLE :: wgfe(:), wgfn(:), wgft(:)
!
      REAL*8, ALLOCATABLE :: usfe(:,:), usfn(:,:), usft(:,:)
      REAL*8, ALLOCATABLE :: vsfe(:,:), vsfn(:,:), vsft(:,:)
      REAL*8, ALLOCATABLE :: wsfe(:,:), wsfn(:,:), wsft(:,:)
!
      REAL*8, DIMENSION(:), ALLOCATABLE  :: kpgv
      REAL*8, DIMENSION(:,:), ALLOCATABLE  :: appu, appv, appw
!
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
      SUBROUTINE fieldn
!----------------------------------------------------------------------
! ... Compute explicitly and store all fields and physical parameters
! ... at time ndt
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
         CALL subscr(ijk)
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

!       END IF
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
      USE convective_fluxes_u, ONLY: flu
      USE convective_fluxes_v, ONLY: flv
      USE convective_fluxes_w, ONLY: flw
      USE control_flags, ONLY: job_type
      USE gas_solid_density, ONLY: rog, rgp, rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE grid, ONLY: dx, dy, dz, fl_l
      USE grid, ONLY: indx, indy, indz, inx, inxb
      USE momentum_transfer, ONLY: kdrags, inter
      USE pressure_epsilon, ONLY: ep, p
      USE time_parameters, ONLY: dt
      USE turbulence_model, ONLY: iss, iturb
      USE gas_solid_viscosity, ONLY: viscg, viscs
      USE gas_solid_viscosity, ONLY: mug
      USE gas_solid_viscosity, ONLY: gvisx, gvisy, gvisz, pvisx, pvisy, pvisz
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm, ijkt
      USE set_indexes, ONLY: stencil, nb, rnb
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
      TYPE(stencil) :: u, v, w, dens
      REAL*8, ALLOCATABLE :: nul(:)
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
      CALL viscg        
!
! ... Calculate particles viscous stress tensor
!
      CALL viscs
!
! ... Allocate and initialize local arrays (gas).
!
      ALLOCATE(rug(ncdom), rwg(ncdom))
      ALLOCATE(ugfe(ncdom), ugft(ncdom))
      ALLOCATE(wgfe(ncdom), wgft(ncdom))
      rug = 0.0D0
      rwg = 0.0D0
      ugfe = 0.0D0; ugft = 0.0D0
      ugfw = 0.0D0; ugfb = 0.0D0
      wgfe = 0.0D0; wgft = 0.0D0
      wgfw = 0.0D0; wgfb = 0.0D0

      IF (job_type == '3D') THEN
        ALLOCATE(rvg(ncdom))
        ALLOCATE(ugfn(ncdom))
        ALLOCATE(vgfe(ncdom), vgfn(ncdom), vgft(ncdom))
        ALLOCATE(wgfn(ncdom))
        rvg = 0.0D0
        ugfn = 0.0D0
        ugfs = 0.0D0
        vgfe = 0.0D0; vgfn = 0.0D0; vgft = 0.0D0
        vgfw = 0.0D0; vgfs = 0.0D0; vgfb = 0.0D0
        wgfn = 0.0D0
        wgfs = 0.0D0
      END IF
!
! ... Allocate and initialize local arrays (particles).
!
      ALLOCATE(rus(ncdom,nsolid), rws(ncdom,nsolid))
      ALLOCATE(usfe(ncdom,nsolid), usft(ncdom,nsolid))
      ALLOCATE(wsfe(ncdom,nsolid), wsft(ncdom,nsolid))

      rus = 0.0D0
      rws = 0.0D0
      usfe = 0.0D0; usft = 0.0D0
      usfw = 0.0D0; usfb = 0.0D0
      wsfe = 0.0D0; wsft = 0.0D0
      wsfw = 0.0D0; wsfb = 0.0D0

      IF (job_type == '3D') THEN
        ALLOCATE(rvs(ncdom,nsolid))
        ALLOCATE(usfn(ncdom,nsolid))
        ALLOCATE(vsfe(ncdom,nsolid), vsfn(ncdom,nsolid), vsft(ncdom,nsolid))
        ALLOCATE(wsfn(ncdom,nsolid))
        rvs = 0.0D0
        usfn = 0.0D0
        usfs = 0.0D0
        vsfe = 0.0D0; vsfn = 0.0D0; vsft = 0.0D0
        vsfw = 0.0D0; vsfs = 0.0D0; vsfb = 0.0D0
        wsfn = 0.0D0
        wsfs = 0.0D0
      END IF
!
! ... Allocate and initialize interphase terms.
!
      ALLOCATE(kpgv(nsolid))
      ALLOCATE(appu(ncdom, ((nsolid+1)**2+(nsolid+1))/2),   &
               appw(ncdom, ((nsolid+1)**2+(nsolid+1))/2))

      kpgv = 0.0D0
      appu = 0.0D0
      appw = 0.0D0
      IF (job_type == '3D') THEN
        ALLOCATE(appv(ncdom, ((nsolid+1)**2+(nsolid+1))/2) )
        appv = 0.0D0
      END IF
      ALLOCATE(nul( ((nsolid+1)**2+(nsolid+1))/2) )
      nul = 0.D0
!
! ... Compute fluxes on East, North and Top sides of a cell
! ... in the whole computational domain.
!
      DO ijk = 1, ncint
        IF(fl_l(ijk) == 1) THEN
          CALL subscr(ijk)
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

            CALL nb(dens,rgp,ijk)
            CALL rnb(u,ug,ijk)
            CALL rnb(v,vg,ijk)
            CALL rnb(w,wg,ijk)

            CALL flu(ugfe(ijk), ugfn(ijk), ugft(ijk),                     &
                    ugfe(imjk), ugfn(ijmk), ugft(ijkm), dens, u, v, w, ijk)
            CALL flv(vgfe(ijk), vgfn(ijk), vgft(ijk),                        &
                    vgfe(imjk), vgfn(ijmk), vgft(ijkm), dens, u, v, w, ijk)
            CALL flw(wgfe(ijk), wgfn(ijk), wgft(ijk),                        &
                    wgfe(imjk), wgfn(ijmk), wgft(ijkm), dens, u, v, w, ijk)
!
            DO is = 1, nsolid
              CALL nb(dens,rlk(:,is),ijk)
              CALL rnb(u,us(:,is),ijk)
              CALL rnb(v,vs(:,is),ijk)
              CALL rnb(w,ws(:,is),ijk)
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
! ... fluxes on West, South and Bottom sides keep values 
! ... of East, North and Top fluxes from neighbouring cells.
!
      DO ijk = 1, ncint
        IF(fl_l(ijk) == 1) THEN
          CALL meshinds(ijk,imesh,i,j,k)
          CALL subscr(ijk)
!
! ... West, South and Bottom fluxes (gas)
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
          rug(ijk) = rugn(ijk) + dt * gvisx(ijk)                     &
     &      - dt * indxp * 2.D0 * ugfx * inxb(i)                     &
     &      - dt * indy(j) * ugfy                                    &
     &      - dt * indz(k) * ugfz   
!
          rwg(ijk) = rwgn(ijk) + dt * gvisz(ijk)                     &
     &      + dt * (dz(k+1)*rgp(ijk)+dz(k)*rgp(ijkt))*indzp * gravz    &
     &      - dt * indx(i) * wgfx * inx(i)                           &
     &      - dt * indy(j) * wgfy                                    &
     &      - dt * indzp * 2.D0 * wgfz
!
          IF (job_type == '3D') THEN

            rvg(ijk) = rvgn(ijk) + dt * gvisy(ijk)                     &
     &        - dt * indx(i) * vgfx                                    &
     &        - dt * indyp * 2.D0 * vgfy                               &
     &        - dt * indz(k) * vgfz    

          END IF

! ... same procedure carried out for particulate phases
!
          DO is = 1, nsolid
!
! ... West, South and Bottom fluxes (particles)
! 
            usfw = usfe(imjk,is)
            usfb = usft(ijkm,is)
            wsfw = wsfe(imjk,is)
            wsfb = wsft(ijkm,is)
!
            usfx = usfe(ijk,is) - usfw
            usfz = usft(ijk,is) - usfb
            wsfx = wsfe(ijk,is) - wsfw
            wsfz = wsft(ijk,is) - wsfb
!
            IF (job_type == '2D') THEN
            
              usfy = 0.D0
              vsfx = 0.D0
              vsfy = 0.D0
              vsfz = 0.D0
              wsfy = 0.D0

            ELSE IF (job_type == '3D') THEN

              usfs = usfn(ijmk,is)
              vsfw = vsfe(imjk,is)
              vsfs = vsfn(ijmk,is)
              vsfb = vsft(ijkm,is)
              wsfs = wsfn(ijmk,is)
!
              usfy = usfn(ijk,is) - usfs
              vsfx = vsfe(ijk,is) - vsfw
              vsfy = vsfn(ijk,is) - vsfs
              vsfz = vsft(ijk,is) - vsfb
              wsfy = wsfn(ijk,is) - wsfs

            END IF
!
! ... compute explicit (tilde) terms in the momentum equation (particles)
! 
              rus(ijk,is) = rusn(ijk,is) + dt*pvisx(ijk,is)              &
     &         - dt*indxp*2.D0* usfx * inxb(i)                            &
     &         - dt*indy(j)* usfy                                         &
     &         - dt*indz(k)* usfz                           
!
              rws(ijk,is) = rwsn(ijk,is) + dt*pvisz(ijk,is)              &
     &         + dt*(rlk(ijk,is)*dz(k+1)+rlk(ijkt,is)*dz(k))*indzp*gravz &
     &         - dt*indx(i)* wsfx * inx(i)                               &
     &         - dt*indy(j)* wsfy                                        &
     &         - dt*indzp*2.D0* wsfz                          
!
            IF (job_type == '3D') THEN

              rvs(ijk,is) = rvsn(ijk,is) + dt*pvisy(ijk,is)              &
     &         - dt*indx(i)* vsfx                                        &
     &         - dt*indyp*2.D0* vsfy                                     &
     &         - dt*indz(k)* vsfz                             

            END IF
!
! ... Compute the gas-particle drag coefficients
!
            dugs = ( (ug(ijk)-us(ijk,is)) + (ug(imjk)-us(imjk,is)) )*0.5D0
            dwgs = ( (wg(ijk)-ws(ijk,is)) + (wg(ijkm)-ws(ijkm,is)) )*0.5D0
            IF (job_type == '2D') THEN
              dvgs = 0.D0
            ELSE IF (job_type == '3D') THEN
              dvgs = ( (vg(ijk)-vs(ijk,is)) + (vg(ijmk)-vs(ijmk,is)) )*0.5D0
            END IF

            CALL kdrags(kpgv(is), dugs, dvgs, dwgs, ep(ijk),         &
                        rgp(ijk), rlk(ijk,is), mug(ijk), is)                  
          END DO 
!
! ... Compute the particle-particle coefficients and the interphase matrix
!
          IF (job_type == '2D') THEN
            CALL inter(appu(ijk,:), nul(:), appw(ijk,:), kpgv(:),    &
     &                 us, us, ws, rlk, ijk)
          ELSE IF (job_type == '2D') THEN
            CALL inter(appu(ijk,:), appv(ijk,:), appw(ijk,:), kpgv(:),    &
     &                 us, vs, ws, rlk, ijk)
          END IF
!
        END IF
      END DO
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
      END MODULE tilde_momentum
!----------------------------------------------------------------------
