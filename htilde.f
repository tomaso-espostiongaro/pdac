!-----------------------------------------------------------------------
      MODULE tilde_energy
!-----------------------------------------------------------------------
      IMPLICIT NONE
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: rhg
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rhs
!
      REAL*8, ALLOCATABLE :: egfe(:), egfn(:), egft(:)
      REAL*8, ALLOCATABLE :: esfe(:,:), esfn(:,:), esft(:,:)
!
      REAL*8, ALLOCATABLE :: hgfe(:), hgfn(:), hgft(:)
      REAL*8, ALLOCATABLE :: hsfe(:,:), hsfn(:,:), hsft(:,:)
!
      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE htilde
!
      USE dimensions
      USE convective_fluxes, ONLY: fsc
      USE diffusive_fluxes, ONLY: hotc
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: sieg, sies, tg, ts
      USE gas_solid_viscosity, ONLY: kapg
      USE grid, ONLY: dx, dy, dz, indx, indy, indz
      USE grid, ONLY: fl_l
      USE grid, ONLY: ncint, ncdom, myijk, data_exchange
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE particles_constants, ONLY: inrl, kap
      USE pressure_epsilon, ONLY: p, pn, ep
      USE set_indexes
      USE time_parameters, ONLY: dt
      USE turbulence, ONLY: kapgt, iturb
      IMPLICIT NONE
!
      REAL*8 :: egfw, egfs, egfb, esfw, esfs, esfb
      REAL*8 :: hrexs, hrexg
      REAL*8 :: flx, prs, upxyz
      REAL*8 :: dxc, dyc, dzc, ugc, vgc, wgc
      INTEGER :: is, m, l, is1
      INTEGER :: i, j, k, ijk, imesh
      TYPE(stencil) :: u, v, w, dens, enth
      TYPE(stencil) :: eps, temp, kappa
!
      ALLOCATE(rhg(ncint))
      ALLOCATE(rhs(ncint,nsolid))
!
      ALLOCATE(egfe(ncdom), egfn(ncdom), egft(ncdom))
      ALLOCATE(esfe(ncdom,nsolid), esfn(ncdom,nsolid), esft(ncdom,nsolid))
!
      ALLOCATE(hgfe(ncdom), hgfn(ncdom), hgft(ncdom))
      ALLOCATE(hsfe(ncdom,nsolid), hsfn(ncdom,nsolid), hsft(ncdom,nsolid))
!
      egfe = 0.0D0;  hgfe = 0.0D0
      egfn = 0.0D0;  hgfn = 0.0D0
      egft = 0.0D0;  hgft = 0.0D0
      esfe = 0.0D0;  hsfe = 0.0D0
      esfn = 0.0D0;  hsfn = 0.0D0
      esft = 0.0D0;  hsft = 0.0D0
!
      rhg = 0.0D0
      rhs = 0.0D0      
!
      IF (iturb >= 1) THEN
        kapgt = kapgt + kapg
      ELSE IF (iturb == 0) THEN
        kapgt = kapg
      END IF
!
      CALL data_exchange(sieg)
      CALL data_exchange(sies)
      CALL data_exchange(tg)
      CALL data_exchange(ts)
      CALL data_exchange(kapgt)
!
      DO ijk = 1, ncint
        imesh = myijk( ip0_jp0_kp0_, ijk)
        IF(fl_l(ijk) == 1) THEN
          CALL subscr(ijk)
!
! ... convective fluxes (gas)
!
          CALL nb(enth,sieg,ijk)
          CALL nb(dens,rgp,ijk)
          CALL rnb(u,ug,ijk)
          CALL rnb(v,vg,ijk)
          CALL rnb(w,wg,ijk)

          CALL fsc(egfe(ijk), egfn(ijk), egft(ijk),    &
                   egfe(imjk), egfn(ijmk), egft(ijkm),    &
                   enth, dens, u, v, w, ijk)
!
! ... diffusive fluxes (gas)
!
          CALL nb(eps,ep,ijk)
          CALL nb(temp,tg,ijk)
          CALL nb(kappa,kapgt,ijk)
          CALL hotc(hgfe(ijk), hgfn(ijk), hgft(ijk),      &
                    hgfe(imjk), hgfn(ijmk), hgft(ijkm),      &
                    eps, temp, kappa, ijk)
!
          egfe(ijk) = egfe(ijk) - hgfe(ijk)
          egfn(ijk) = egfn(ijk) - hgfn(ijk)
          egft(ijk) = egft(ijk) - hgft(ijk)
!
          IF (fl_l(imjk) /= 1) egfe(imjk) = egfe(imjk) - hgfe(imjk)
          IF (fl_l(ijmk) /= 1) egfn(ijmk) = egfn(ijmk) - hgfn(ijmk)
          IF (fl_l(ijkm) /= 1) egft(ijkm) = egft(ijkm) - hgft(ijkm)
!
          DO is=1, nsolid
!
! ... convective fluxes (particles)
!
            CALL nb(enth,sies(:,is),ijk)
            CALL nb(dens,rlk(:,is),ijk)
            CALL rnb(u, us(:,is),ijk)
            CALL rnb(v, vs(:,is),ijk)
            CALL rnb(w, ws(:,is),ijk)
            
            CALL fsc(esfe(ijk, is), esfn(ijk, is), esft(ijk, is),  &
                     esfe(imjk, is), esfn(ijmk, is), esft(ijkm, is),  &
                     enth, dens, u, v, w, ijk)
!
! ... diffusive fluxes (particles)
!
            eps   = inrl(is) * dens
            CALL nb(temp,ts(:,is),ijk)
            kappa = cte(kap(is))
            CALL hotc(hsfe(ijk,is), hsfn(ijk, is), hsft(ijk, is),    &
                      hsfe(imjk,is), hsfn(ijmk, is), hsft(ijkm, is),    &
                      eps, temp, kappa, ijk)
!
            esfe(ijk, is) = esfe(ijk, is) - hsfe(ijk,is)
            esfn(ijk, is) = esfn(ijk, is) - hsfn(ijk, is)
            esft(ijk, is) = esft(ijk, is) - hsft(ijk, is)
!
            IF (fl_l(imjk) /= 1) esfe(imjk,is) = esfe(imjk,is) - hsfe(imjk,is)
            IF (fl_l(ijmk) /= 1) esfn(ijmk,is) = esfn(ijmk,is) - hsfn(ijmk,is)
            IF (fl_l(ijkm) /= 1) esft(ijkm,is) = esft(ijkm,is) - hsft(ijkm,is)
!
          END DO
        END IF
      END DO
!
      CALL data_exchange(egfe)
      CALL data_exchange(egfn)
      CALL data_exchange(egft)
      CALL data_exchange(esfe)
      CALL data_exchange(esfn)
      CALL data_exchange(esft)
!
! ... fluxes on left and bottom sides keep values
! ... entering from neighbouring cells.
!
      DO ijk = 1, ncint
        imesh = myijk( ip0_jp0_kp0_, ijk)
        IF(fl_l(ijk) == 1) THEN
          CALL subscr(ijk)
          i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
          j = MOD( imesh - 1, nx*ny ) / nx + 1
          k = ( imesh - 1 ) / ( nx*ny ) + 1
! 
          egfw = egfe(imjk)
          egfs = egfn(ijmk)
          egfb = egft(ijkm)
!
          flx = dt * indx(i) * (egfe(ijk)-egfw)   +   &
                dt * indy(j) * (egfn(ijk)-egfs)   +   &
                dt * indz(k) * (egft(ijk)-egfb)
!
          dxc = dx(i)+(dx(i+1)+dx(i-1))/2.D0 
          dyc = dy(j)+(dy(j+1)+dy(j-1))/2.D0
          dzc = dz(k)+(dz(k+1)+dz(k-1))/2.D0
          ugc = (ug(ijk)+ug(imjk))/2.D0
          vgc = (vg(ijk)+vg(ijmk))/2.D0
          wgc = (wg(ijk)+wg(ijkm))/2.D0
          upxyz= dt/dxc * ugc * (p(ijke)-p(ijkw)) +   &
                 dt/dyc * vgc * (p(ijkn)-p(ijks)) +   &
                 dt/dzc * wgc * (p(ijkt)-p(ijkb))
!
          prs = ep(ijk) * (p(ijk) - pn(ijk) + upxyz)
!
          rhg(ijk) = - flx + prs
!
! ... Same procedure carried out for solids
!
          DO is=1, nsolid
           IF (rlk(imjk,is) * inrl(is) <= 1.D-9) THEN
            esfw = 0.0D0
           ELSE
            esfw = esfe(imjk, is)
           END IF
           IF (rlk(ijmk,is) * inrl(is) <= 1.D-9) THEN
            esfs = 0.0D0
           ELSE
            esfs = esfn(ijmk, is)
           END IF
           IF (rlk(ijkm,is) * inrl(is) <= 1.D-9) THEN
            esfb = 0.0D0
           ELSE
            esfb = esft(ijkm, is)
           END IF
!
            flx = dt * indx(i) * (esfe(ijk, is) - esfw)  +  &
                  dt * indy(j) * (esfn(ijk, is) - esfs)  +  &
                  dt * indz(k) * (esft(ijk, is) - esfb)
!
            rhs(ijk,is) = - flx
          END DO
!
        END IF
      END DO

!
      DEALLOCATE(egfe, egfn, egft)
      DEALLOCATE(esfe, esfn, esft)
!
      DEALLOCATE(hgfe, hgfn, hgft)
      DEALLOCATE(hsfe, hsfn, hsft)
!
      RETURN
      END SUBROUTINE htilde
!-----------------------------------------------------------
      END MODULE tilde_energy
!-----------------------------------------------------------
