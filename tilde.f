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
      SUBROUTINE fieldn
!----------------------------------------------------------------------
! ... Compute explicitly and store all fields at time ndt
!
      USE dimensions
      USE eos_gas, ONLY: rgpgc, rgpgcn, xgc
      USE gas_solid_density, ONLY: rgp, rgpn, rlk, rlkn
      USE gas_solid_temperature, ONLY: sieg, siegn, sies, siesn, tg
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_viscosity, ONLY: viscon, mug, kapg
      USE grid, ONLY: dz, dy, dx, fl_l
      USE grid, ONLY: indx, indy, indz
      USE grid, ONLY: ncint, myijk, ncdom, data_exchange
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
      IF (ALLOCATED(rugn)) DEALLOCATE(rugn)
      IF (ALLOCATED(rvgn)) DEALLOCATE(rvgn)
      IF (ALLOCATED(rwgn)) DEALLOCATE(rwgn)
      IF (ALLOCATED(rusn)) DEALLOCATE(rusn)
      IF (ALLOCATED(rvsn)) DEALLOCATE(rvsn)
      IF (ALLOCATED(rwsn)) DEALLOCATE(rwsn)
      ALLOCATE( rugn(ncint),  rvgn(ncint), rwgn(ncint))
      ALLOCATE( rusn(nsolid,ncint), rvsn(nsolid,ncint), rwsn(nsolid,ncint))
      rugn = 0.D0; rvgn = 0.D0; rwgn = 0.D0
      rusn = 0.D0; rvsn = 0.D0; rwsn = 0.D0
!
      CALL data_exchange(rgp)
      CALL data_exchange(rlk)
!
      DO ijk = 1, ncint
       imesh = myijk( ip0_jp0_kp0_,ijk)
       IF(fl_l(ijk) == 1) THEN
         CALL subscr(ijk)
         i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
         j = MOD( imesh - 1, nx*ny ) / nx + 1
         k = ( imesh - 1 ) / ( nx*ny ) + 1
!
         dxp=dx(i)+dx(i+1)
         dyp=dy(j)+dy(j+1)
         dzp=dz(k)+dz(k+1)
         indxp=1.D0/dxp
         indyp=1.D0/dyp
         indzp=1.D0/dzp
!
         rgp_e = (dx(i+1)*rgp(ijk)+dx(i)*rgp(ijke))*indxp
         rgp_n = (dy(j+1)*rgp(ijk)+dy(j)*rgp(ijkn))*indyp
         rgp_t = (dz(k+1)*rgp(ijk)+dz(k)*rgp(ijkt))*indzp
!
         rugn(ijk)  = rgp_e * ug(ijk)
         rvgn(ijk)  = rgp_n * vg(ijk)
         rwgn(ijk)  = rgp_t * wg(ijk)
!
         pn(ijk)    = p(ijk)
         rgpn(ijk)  = rgp(ijk)
         siegn(ijk) = sieg(ijk)
         DO ig=1,ngas
           rgpgcn(ig,ijk) = rgpgc(ig,ijk)
         END DO
!
         DO is = 1, nsolid
          rlk_e = (rlk(is,ijk)*dx(i+1)+rlk(is,ijke)*dx(i))*indxp
          rlk_n = (rlk(is,ijk)*dy(j+1)+rlk(is,ijkn)*dy(j))*indyp
          rlk_t = (rlk(is,ijk)*dz(k+1)+rlk(is,ijkt)*dz(k))*indzp
!
          rusn(is,ijk)  = rlk_e * us(is,ijk)
          rvsn(is,ijk)  = rlk_n * vs(is,ijk)
          rwsn(is,ijk)  = rlk_t * ws(is,ijk)
!
          rlkn(is,ijk)  = rlk(is,ijk)
          siesn(is,ijk) = sies(is,ijk)
         END DO
!
! ... Compute the temperature-dependent gas viscosity and th. conductivity
!
         CALL viscon(mug(ijk), kapg(ijk), xgc(:,ijk), tg(ijk))
!	
       END IF
      END DO
!
      RETURN
      END SUBROUTINE fieldn
!----------------------------------------------------------------------
      SUBROUTINE tilde
!
      USE atmosphere, ONLY: gravz
      USE grid, ONLY: fl_l
      USE dimensions
      USE convective_fluxes, ONLY: flu, flv, flw
      USE gas_solid_density, ONLY: rog, rgp, rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE grid, ONLY: dx, dy, dz
      USE grid, ONLY: indx, indy, indz
      USE momentum_transfer, ONLY: kdrags, inter
      USE pressure_epsilon, ONLY: ep, p
      USE time_parameters, ONLY: dt
      USE turbulence, ONLY: iss, iturb
      USE gas_solid_viscosity, ONLY: viscg, viscs
      USE gas_solid_viscosity, ONLY: mug
      USE gas_solid_viscosity, ONLY: gvisx, gvisy, gvisz, pvisx, pvisy, pvisz
      USE grid, ONLY: ncint, myijk, ncdom, data_exchange
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes
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
!
      ALLOCATE(gvisx(ncint), gvisy(ncint), gvisz(ncint))
      ALLOCATE(pvisx(nsolid, ncint), pvisy(nsolid,ncint), pvisz(nsolid, ncint))
      gvisx = 0.D0; gvisy = 0.D0; gvisz = 0.D0
      pvisx = 0.D0; pvisy = 0.D0; pvisz = 0.D0

      IF (iturb == 0) THEN
        CALL data_exchange(ug)
        CALL data_exchange(vg)
        CALL data_exchange(wg)
      END IF
      IF (iss == 0) THEN
        CALL data_exchange(us)
        CALL data_exchange(vs)
        CALL data_exchange(ws)
      END IF
      CALL data_exchange(ep)
!
! ... Calculate gvisx, gvisy and gvisz (gas viscous stress tensor).
!
      CALL viscg        
!
! ... Calculate pvisx, pvisy and pvisz (particles viscous stress tensor).
!
      CALL viscs
!
! ... Allocate and initialize local arrays (gas).
!
      ALLOCATE(rug(ncdom),  rvg(ncdom),  rwg(ncdom))
      ALLOCATE(ugfe(ncdom), ugfn(ncdom), ugft(ncdom))
      ALLOCATE(vgfe(ncdom), vgfn(ncdom), vgft(ncdom))
      ALLOCATE(wgfe(ncdom), wgfn(ncdom), wgft(ncdom))

      rug = 0.0D0
      rvg = 0.0D0
      rwg = 0.0D0
      ugfe = 0.0D0; ugfn = 0.0D0; ugft = 0.0D0
      ugfw = 0.0D0; ugfs = 0.0D0; ugfb = 0.0D0
      vgfe = 0.0D0; vgfn = 0.0D0; vgft = 0.0D0
      vgfw = 0.0D0; vgfs = 0.0D0; vgfb = 0.0D0
      wgfe = 0.0D0; wgfn = 0.0D0; wgft = 0.0D0
      wgfw = 0.0D0; wgfs = 0.0D0; wgfb = 0.0D0
!
! ... Allocate and initialize local arrays (particles).
!
      ALLOCATE(rus(nsolid,ncdom),  rvs(nsolid,ncdom),  rws(nsolid,ncdom))
      ALLOCATE(usfe(nsolid,ncdom), usfn(nsolid,ncdom), usft(nsolid,ncdom))
      ALLOCATE(vsfe(nsolid,ncdom), vsfn(nsolid,ncdom), vsft(nsolid,ncdom))
      ALLOCATE(wsfe(nsolid,ncdom), wsfn(nsolid,ncdom), wsft(nsolid,ncdom))

      rus = 0.0D0
      rvs = 0.0D0
      rws = 0.0D0
      usfe = 0.0D0; usfn = 0.0D0; usft = 0.0D0
      usfw = 0.0D0; usfs = 0.0D0; usfb = 0.0D0
      vsfe = 0.0D0; vsfn = 0.0D0; vsft = 0.0D0
      vsfw = 0.0D0; vsfs = 0.0D0; vsfb = 0.0D0
      wsfe = 0.0D0; wsfn = 0.0D0; wsft = 0.0D0
      wsfw = 0.0D0; wsfs = 0.0D0; wsfb = 0.0D0
!
! ... Allocate and initialize interphase terms.
!
      ALLOCATE(kpgv(nsolid))
      ALLOCATE(appu(((nsolid+1)**2+(nsolid+1))/2, ncdom),   &
               appv(((nsolid+1)**2+(nsolid+1))/2, ncdom),   &
               appw(((nsolid+1)**2+(nsolid+1))/2, ncdom))

      kpgv = 0.0D0
      appu = 0.0D0
      appv = 0.0D0
      appw = 0.0D0
!
! ... Compute fluxes on East, North and Top sides of a cell
! ... in the whole computational domain.
!
      DO ijk = 1, ncint
        imesh = myijk( ip0_jp0_kp0_, ijk)
        IF(fl_l(ijk) == 1) THEN
          CALL subscr(ijk)
          dens = nb(rgp,ijk)
          u    = rnb(ug,ijk)
          v    = rnb(vg,ijk)
          w    = rnb(wg,ijk)
!
          CALL flu(ugfe(ijk), ugfn(ijk), ugft(ijk),                        &
                  ugfe(imjk), ugfn(ijmk), ugft(ijkm), dens, u, v, w, ijk)
          CALL flv(vgfe(ijk), vgfn(ijk), vgft(ijk),                        &
                  vgfe(imjk), vgfn(ijmk), vgft(ijkm), dens, u, v, w, ijk)
          CALL flw(wgfe(ijk), wgfn(ijk), wgft(ijk),                        &
                  wgfe(imjk), wgfn(ijmk), wgft(ijkm), dens, u, v, w, ijk)
!
          DO is = 1, nsolid
            dens = nb(rlk,is,ijk)
            u    = rnb(us,is,ijk)
            v    = rnb(vs,is,ijk)
            w    = rnb(ws,is,ijk)
!
            CALL flu(usfe(is,ijk), usfn(is,ijk), usft(is,ijk),           &
                    usfe(is,imjk), usfn(is,ijmk), usft(is,ijkm),        &
                    dens, u, v, w, ijk)
            CALL flv(vsfe(is,ijk), vsfn(is,ijk), vsft(is,ijk),           &
                    vsfe(is,imjk), vsfn(is,ijmk), vsft(is,ijkm),        &
                    dens, u, v, w, ijk)
            CALL flw(wsfe(is,ijk), wsfn(is,ijk), wsft(is,ijk),           &
                    wsfe(is,imjk), wsfn(is,ijmk), wsft(is,ijkm),        &
                    dens, u, v, w, ijk)
!
          END DO
        END IF         
      END DO
!
      CALL data_exchange(ugfe)
      CALL data_exchange(ugfn)
      CALL data_exchange(ugft)
      CALL data_exchange(vgfe)
      CALL data_exchange(vgfn)
      CALL data_exchange(vgft)
      CALL data_exchange(wgfe)
      CALL data_exchange(wgfn)
      CALL data_exchange(wgft)

      CALL data_exchange(usfe)
      CALL data_exchange(usfn)
      CALL data_exchange(usft)
      CALL data_exchange(vsfe)
      CALL data_exchange(vsfn)
      CALL data_exchange(vsft)
      CALL data_exchange(wsfe)
      CALL data_exchange(wsfn)
      CALL data_exchange(wsft)
!
! ... fluxes on West, South and Bottom sides keep values 
! ... of East, North and Top fluxes from neighbouring cells.
!
      DO ijk = 1, ncint
        imesh = myijk( ip0_jp0_kp0_, ijk)
        IF(fl_l(ijk) == 1) THEN
          i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
          j = MOD( imesh - 1, nx*ny ) / nx + 1
          k = ( imesh - 1 ) / ( nx*ny ) + 1
          CALL subscr(ijk)
!
          dxp=dx(i)+dx(i+1)
          dyp=dy(j)+dy(j+1)
          dzp=dz(k)+dz(k+1)
          indxp=1.D0/dxp
          indyp=1.D0/dyp
          indzp=1.D0/dzp
!
! ... West, South and Bottom fluxes (gas)
!
          ugfw = ugfe(imjk)
          ugfs = ugfn(ijmk)
          ugfb = ugft(ijkm)
          vgfw = vgfe(imjk)
          vgfs = vgfn(ijmk)
          vgfb = vgft(ijkm)
          wgfw = wgfe(imjk)
          wgfs = wgfn(ijmk)
          wgfb = wgft(ijkm)
!
! ... compute explicit (tilde) terms in the momentum equation (gas)
! 
          ugfx = ugfe(ijk) - ugfw
          ugfy = ugfn(ijk) - ugfs
          ugfz = ugft(ijk) - ugfb
!
          rug(ijk) = rugn(ijk) + dt * gvisx(ijk)                     &
     &      - dt * indxp * 2.D0 * ugfx                               &
     &      - dt * indy(j) * ugfy                                    &
     &      - dt * indz(k) * ugfz   
!
          vgfx = vgfe(ijk) - vgfw
          vgfy = vgfn(ijk) - vgfs
          vgfz = vgft(ijk) - vgfb
!
          rvg(ijk) = rvgn(ijk) + dt * gvisy(ijk)                     &
     &      - dt * indx(i) * vgfx                                    &
     &      - dt * indyp * 2.D0 * vgfy                               &
     &      - dt * indz(k) * vgfz    
!
          wgfx = wgfe(ijk) - wgfw
          wgfy = wgfn(ijk) - wgfs
          wgfz = wgft(ijk) - wgfb
!
          rwg(ijk) = rwgn(ijk) + dt * gvisz(ijk)                     &
     &      + dt * (dz(k+1)*rgp(ijk)+dz(k)*rgp(ijkt))*indzp * gravz    &
     &      - dt * indx(i) * wgfx                                    &
     &      - dt * indy(j) * wgfy                                    &
     &      - dt * indzp * 2.D0 * wgfz
!
! ... same procedure carried out for particulate phases
!
          DO is = 1, nsolid
!
! ... West, South and Bottom fluxes (particles)
! 
            usfw = usfe(is,imjk)
            usfs = usfn(is,ijmk)
            usfb = usft(is,ijkm)
            vsfw = vsfe(is,imjk)
            vsfs = vsfn(is,ijmk)
            vsfb = vsft(is,ijkm)
            wsfw = wsfe(is,imjk)
            wsfs = wsfn(is,ijmk)
            wsfb = wsft(is,ijkm)
!
! ... compute explicit (tilde) terms in the momentum equation (particles)
! 
              usfx = usfe(is,ijk) - usfw
              usfy = usfn(is,ijk) - usfs
              usfz = usft(is,ijk) - usfb
!
              rus(is,ijk) = rusn(is,ijk) + dt*pvisx(is,ijk)              &
     &         - dt*indxp*2.D0* usfx                                      &
     &         - dt*indy(j)* usfy                                         &
     &         - dt*indz(k)* usfz                           
!
              vsfx = vsfe(is,ijk) - vsfw
              vsfy = vsfn(is,ijk) - vsfs
              vsfz = vsft(is,ijk) - vsfb
!
              rvs(is,ijk) = rvsn(is,ijk) + dt*pvisy(is,ijk)              &
     &         - dt*indx(i)* vsfx                                        &
     &         - dt*indyp*2.D0* vsfy                                     &
     &         - dt*indz(k)* vsfz                             
!
              wsfx = wsfe(is,ijk) - wsfw
              wsfy = wsfn(is,ijk) - wsfs
              wsfz = wsft(is,ijk) - wsfb
!
              rws(is,ijk) = rwsn(is,ijk) + dt*pvisz(is,ijk)              &
     &         + dt*(rlk(is,ijk)*dz(k+1)+rlk(is,ijkt)*dz(k))*indzp*gravz &
     &         - dt*indx(i)* wsfx                                        &
     &         - dt*indy(j)* wsfy                                        &
     &         - dt*indzp*2.D0* wsfz                          
!
! ... Compute the gas-particle drag coefficients
!
            dugs = ( (ug(ijk)-us(is,ijk)) + (ug(imjk)-us(is,imjk)) )*0.5D0
            dvgs = ( (vg(ijk)-vs(is,ijk)) + (vg(ijmk)-vs(is,ijmk)) )*0.5D0
            dwgs = ( (wg(ijk)-ws(is,ijk)) + (wg(ijkm)-ws(is,ijkm)) )*0.5D0

            CALL kdrags(kpgv(is), dugs, dvgs, dwgs, ep(ijk),         &
                        rgp(ijk), rlk(is,ijk), mug(ijk), is)                  

          END DO 
!
! ... Compute the particle-particle coefficients and the interphase matrix
!
          CALL inter(appu(:,ijk), appv(:,ijk), appw(:,ijk), kpgv(:),    &
     &               us, vs, ws, rlk, ijk)
!
        END IF
      END DO
!      appu = 0.D0
!      appv = 0.D0
!      appw = 0.D0
!
      DEALLOCATE(ugfe, ugfn, ugft)
      DEALLOCATE(vgfe, vgfn, vgft)
      DEALLOCATE(wgfe, wgfn, wgft)
      DEALLOCATE(gvisx, gvisy, gvisz)
!
      DEALLOCATE(usfe, usfn, usft)
      DEALLOCATE(vsfe, vsfn, vsft)
      DEALLOCATE(wsfe, wsfn, wsft)
      DEALLOCATE(pvisx, pvisy, pvisz)
!
      DEALLOCATE(kpgv)
!
      CALL data_exchange(appu)
      CALL data_exchange(appv)
      CALL data_exchange(appw)
      CALL data_exchange(rug)
      CALL data_exchange(rvg)
      CALL data_exchange(rwg)
      CALL data_exchange(rus)
      CALL data_exchange(rvs)
      CALL data_exchange(rws)
!
      RETURN
      END SUBROUTINE tilde
!----------------------------------------------------------------------
      END MODULE tilde_momentum
!----------------------------------------------------------------------
