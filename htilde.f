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
      ALLOCATE(rhs(nsolid,ncint))
!
      ALLOCATE(egfe(ncdom), egfn(ncdom), egft(ncdom))
      ALLOCATE(esfe(nsolid,ncdom), esfn(nsolid,ncdom), esft(nsolid,ncdom))
!
      ALLOCATE(hgfe(ncdom), hgfn(ncdom), hgft(ncdom))
      ALLOCATE(hsfe(nsolid,ncdom), hsfn(nsolid,ncdom), hsft(nsolid,ncdom))
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
      ELSE
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
          enth  = nb(sieg,ijk)
          dens  = nb(rgp,ijk)
          u     = rnb(ug,ijk)
          v     = rnb(vg,ijk)
          w     = rnb(wg,ijk)
          CALL fsc(egfe(ijk), egfn(ijk), egft(ijk),    &
                   egfe(imjk), egfn(ijmk), egft(ijkm),    &
                   enth, dens, u, v, w, ijk)
!
! ... diffusive fluxes (gas)
!
          eps   = nb(ep,ijk)
          temp  = nb(tg,ijk)
          kappa = nb(kapgt,ijk)
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
            enth  = nb(sies,is,ijk)
            dens  = nb(rlk,is,ijk)
            u     = rnb(us,is,ijk)
            v     = rnb(vs,is,ijk)
            w     = rnb(ws,is,ijk)
            CALL fsc(esfe(is, ijk), esfn(is, ijk), esft(is, ijk),  &
                     esfe(is, imjk), esfn(is, ijmk), esft(is, ijkm),  &
                     enth, dens, u, v, w, ijk)
!
! ... diffusive fluxes (particles)
!
            eps   = inrl(is) * nb(rlk,is,ijk)
            temp  = nb(ts,is,ijk)
            kappa = cte(kap(is))
            CALL hotc(hsfe(is,ijk), hsfn(is, ijk), hsft(is,ijk),    &
                      hsfe(is,imjk), hsfn(is, ijmk), hsft(is,ijkm),    &
                      eps, temp, kappa, ijk)
!
            esfe(is, ijk) = esfe(is, ijk) - hsfe(is, ijk)
            esfn(is, ijk) = esfn(is, ijk) - hsfn(is, ijk)
            esft(is, ijk) = esft(is, ijk) - hsft(is, ijk)
!
            IF (fl_l(imjk) /= 1) esfe(is,imjk) = esfe(is,imjk) - hsfe(is,imjk)
            IF (fl_l(ijmk) /= 1) esfn(is,ijmk) = esfn(is,ijmk) - hsfn(is,ijmk)
            IF (fl_l(ijkm) /= 1) esft(is,ijkm) = esft(is,ijkm) - hsft(is,ijkm)
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
           IF (rlk(is,imjk) * inrl(is) <= 1.D-9) THEN
            esfw = 0.0D0
           ELSE
            esfw = esfe(is, imjk)
           END IF
           IF (rlk(is,ijmk) * inrl(is) <= 1.D-9) THEN
            esfs = 0.0D0
           ELSE
            esfs = esfn(is, ijmk)
           END IF
           IF (rlk(is,ijkm) * inrl(is) <= 1.D-9) THEN
            esfb = 0.0D0
           ELSE
            esfb = esft(is, ijkm)
           END IF
!
            flx = dt * indx(i) * (esfe(is,ijk) - esfw)  +  &
                  dt * indy(j) * (esfn(is,ijk) - esfs)  +  &
                  dt * indz(k) * (esft(is,ijk) - esfb)
!
            rhs(is, ijk) = - flx
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
