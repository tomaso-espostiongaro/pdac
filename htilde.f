!-----------------------------------------------------------------------
      MODULE tilde_energy
!-----------------------------------------------------------------------
      IMPLICIT NONE
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: rhg
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rhk
!
      REAL*8, ALLOCATABLE :: egfr(:), egft(:)
      REAL*8, ALLOCATABLE :: egfl(:), egfb(:)
      REAL*8, ALLOCATABLE :: elfr(:,:), elft(:,:)
      REAL*8, ALLOCATABLE :: elfl(:,:), elfb(:,:)
!
      REAL*8, ALLOCATABLE :: hfgr(:), hfgt(:)
      REAL*8, ALLOCATABLE :: hfgl(:), hfgb(:)
      REAL*8, ALLOCATABLE :: hflr(:,:), hflt(:,:)
      REAL*8, ALLOCATABLE :: hfll(:,:), hflb(:,:)
!
      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE htilde
!
      USE dimensions
      USE eulerian_flux, ONLY: fsc_rt, fsc_lb
      USE gas_solid_velocity, ONLY: ug, wg, us, ws
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: sieg, sies, tg, ts
      USE gas_solid_viscosity, ONLY: kapg
      USE grid, ONLY: r, rb, dr, zb, dz, inr, indr, indz 
      USE grid, ONLY: fl_l
      USE grid, ONLY: ncint, ncdom, myijk, data_exchange
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE heat_diffusion, ONLY: hotcg, hotck
      USE particles_constants, ONLY: inrl
      USE pressure_epsilon, ONLY: p, pn, ep
      USE time_parameters, ONLY: dt
      USE set_indexes
      USE turbulence, ONLY: kapgt, iturb
      IMPLICIT NONE
!
      REAL*8 :: c3, hrexs, hrexg, c2, upxy, c1
      REAL*8 :: drc, dzc, ugb, wgb
      INTEGER :: is, m, l, is1, i, j, imesh
      INTEGER :: ij
      TYPE(stencil) :: u, w, dens, enth
      TYPE(stencil) :: eps, temp, kap
!
      ALLOCATE(rhg(ncint))
      ALLOCATE(rhk(nsolid,ncint))
!
      ALLOCATE(egfr(ncdom), egft(ncdom))
      ALLOCATE(egfl(ncdom), egfb(ncdom))
      ALLOCATE(elfr(nsolid,ncdom), elft(nsolid,ncdom))
      ALLOCATE(elfl(nsolid,ncdom), elfb(nsolid,ncdom))
!
      ALLOCATE(hfgr(ncdom), hfgt(ncdom))
      ALLOCATE(hfgl(ncdom), hfgb(ncdom))
      ALLOCATE(hflr(nsolid,ncdom), hflt(nsolid,ncdom))
      ALLOCATE(hfll(nsolid,ncdom), hflb(nsolid,ncdom))
!
      egfr = 0.0D0;  hfgr = 0.0D0
      egft = 0.0D0;  hfgt = 0.0D0
      egfl = 0.0D0;  hfgl = 0.0D0
      egfb = 0.0D0;  hfgb = 0.0D0
      elfr = 0.0D0;  hflr = 0.0D0
      elft = 0.0D0;  hflt = 0.0D0
      elfl = 0.0D0;  hfll = 0.0D0
      elfb = 0.0D0;  hflb = 0.0D0
!
      rhg = 0.0D0
      rhk = 0.0D0      
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
      DO ij = 1, ncint
        imesh = myijk( ip0_jp0_kp0_, ij)
        IF(fl_l(ij).EQ.1) THEN
          CALL subscr(ij)
!
! ... convective radial fluxes (gas)
!
          enth  = nb(sieg,ij)
          dens  = nb(rgp,ij)
          u     = rnb(ug,ij)
          w     = rnb(wg,ij)
          CALL fsc_rt(egfr(ij), egft(ij), enth, dens, u, w, ij)
!
! ... diffusive radial fluxes (gas)
!
          eps  = nb(ep,ij)
          temp = nb(tg,ij)
          kap  = nb(kapgt,ij)
          CALL hotcg(hfgr(ij), hfgt(ij), hfgl(ij), hfgb(ij),  &
               eps, temp, kap, ij)
!
          egfr(ij) = egfr(ij) - hfgr(ij)
          egft(ij) = egft(ij) - hfgt(ij)
!
          DO is=1, nsolid
!
! ... convective radial fluxes (particles)
!
            enth  = nb(sies,is,ij)
            dens  = nb(rlk,is,ij)
            u     = rnb(us,is,ij)
            w     = rnb(ws,is,ij)
            CALL fsc_rt(elfr(is, ij), elft(is, ij), enth, dens, u, w, ij)
!
! ... diffusive radial fluxes (particles)
!
            dens = nb(rlk,is,ij)
            temp = nb(ts,is,ij)
            CALL hotck(hflr(is,ij), hflt(is,ij), hfll(is,ij), hflb(is,ij),  &
                 dens, temp, ij, is)
!
          elfr(is, ij) = elfr(is, ij) - hflr(is, ij)
          elft(is, ij) = elft(is, ij) - hflt(is, ij)
!
          END DO
        END IF
      END DO
!
      CALL data_exchange(egfr)
      CALL data_exchange(egft)
      CALL data_exchange(elfr)
      CALL data_exchange(elft)
!
! ... fluxes on left and bottom sides keep values
! ... entering from neighbouring cells.
! ... on boundaries, fluxes on left and bottom sides must be calculated
!
      DO ij = 1, ncint
        imesh = myijk( ip0_jp0_kp0_, ij)
        IF(fl_l(ij).EQ.1) THEN
          CALL subscr(ij)
          j = ( imesh - 1 ) / nr + 1
          i = MOD( ( imesh - 1 ), nr) + 1
! 
          egfl(ij) = egfr(imj)
          egfb(ij) = egft(ijm)
!
! ... convective vertical fluxes (gas)
!
          enth  = nb(sieg,ij)
          dens  = nb(rgp,ij)
          u     = rnb(ug,ij)
          w     = rnb(wg,ij)
          CALL fsc_lb(egfl(ij), egfb(ij), enth, dens, u, w, ij)
!
! ... diffusive vertical fluxes (gas) (WHY COMPUTED AGAIN???)
!
          eps  = nb(ep,ij)
          temp = nb(tg,ij)
          kap  = nb(kapgt,ij)
          CALL hotcg(hfgr(ij), hfgt(ij), hfgl(ij), hfgb(ij),  &
               eps, temp, kap, ij)
!
            egfl(ij) = egfl(ij) - hfgl(ij)
            egfb(ij) = egfb(ij) - hfgb(ij)
!
          DO is=1, nsolid
           IF (rlk(is,imj) * inrl(is) .LE. 1.D-9) THEN
            elfl(is, ij) = 0.0D0
           ELSE
            elfl(is, ij) = elfr(is, imj)
           END IF
           IF (rlk(is,ijm) * inrl(is) .LE. 1.D-9) THEN
            elfb(is, ij) = 0.0D0
           ELSE
            elfb(is, ij) = elft(is, ijm)
           END IF
!
! ... convective vertical fluxes (particles)
!
            enth  = nb(sies,is,ij)
            dens  = nb(rlk,is,ij)
            u     = rnb(us,is,ij)
            w     = rnb(ws,is,ij)
            CALL fsc_lb(elfl(is, ij), elfb(is, ij), enth, dens, u, w, ij)
!
! ... diffusive vertical fluxes (particles) (WHY COMPUTED AGAIN???)
!
            dens = nb(rlk,is,ij)
            temp = nb(ts,is,ij)
            CALL hotck(hflr(is,ij), hflt(is,ij), hfll(is,ij), hflb(is,ij),  &
                 dens, temp, ij, is)
!
          elfl(is, ij) = elfl(is, ij) - hfll(is, ij)
          elfb(is, ij) = elfb(is, ij) - hflb(is, ij)
!
          END DO
!
          c1 = dt * inr(i)*indr(i) * (egfr(ij)-egfl(ij))     &
             + dt * indz(j)        * (egft(ij)-egfb(ij))
!
          drc = dr(i)+(dr(i+1)+dr(i-1))/2.D0 
          dzc = dz(j)+(dz(j+1)+dz(j-1))/2.D0
          ugb = (ug(ij)+ug(imj))/2.D0
          wgb = (wg(ij)+wg(ijm))/2.D0
          upxy= dt/drc * ugb * (p(ijr)-p(ijl)) + dt/dzc * wgb * (p(ijt)-p(ijb))
!
          c2 = ep(ij) * (p(ij)-pn(ij)+upxy)
!
          rhg(ij) = - c1 + c2
!
          DO is=1,nsolid
            c3 = dt * inr(i)*indr(i) * (elfr(is,ij) - elfl(is,ij))    &
               + dt * indz(j)        * (elft(is,ij) - elfb(is,ij))
!
            rhk(is, ij) = - c3
          END DO
!
        END IF
      END DO
!
      DEALLOCATE(egfr, egft)
      DEALLOCATE(egfl, egfb)
      DEALLOCATE(elfr, elft)
      DEALLOCATE(elfl, elfb)
!
      DEALLOCATE(hfgr, hfgt)
      DEALLOCATE(hfgl, hfgb)
      DEALLOCATE(hflr, hflt)
      DEALLOCATE(hfll, hflb)
!
      RETURN
      END SUBROUTINE htilde
!-----------------------------------------------------------
      END MODULE tilde_energy
