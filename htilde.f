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
      USE gas_solid_velocity, ONLY: ug, vg, uk, vk
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: sieg, siek, tg, tk
      USE gas_solid_viscosity, ONLY: kapg
      USE grid, ONLY: r, rb, dr, zb, dz, inr, indr, indz 
      USE grid, ONLY: fl_l
      USE grid, ONLY: nij_l, nijx_l, myij, data_exchange
      USE heat_diffusion, ONLY: hotcg, hotck
      USE particles_constants, ONLY: inrl
      USE pressure_epsilon, ONLY: p, pn, ep
      USE time_parameters, ONLY: dt
      USE set_indexes
      USE turbulence, ONLY: kapgt, iturb
      IMPLICIT NONE
!
      REAL*8 :: c3, hrexs, hrexg, c2, upxy, c1
      REAL*8 :: drc, dzc, ugb, vgb
      INTEGER :: k, m, l, k1, i, j, imesh
      INTEGER :: ij
!
      ALLOCATE(rhg(nij_l))
      ALLOCATE(rhk(nsolid,nij_l))
!
      ALLOCATE(egfr(nijx_l), egft(nijx_l))
      ALLOCATE(egfl(nijx_l), egfb(nijx_l))
      ALLOCATE(elfr(nsolid,nijx_l), elft(nsolid,nijx_l))
      ALLOCATE(elfl(nsolid,nijx_l), elfb(nsolid,nijx_l))
!
      ALLOCATE(hfgr(nijx_l), hfgt(nijx_l))
      ALLOCATE(hfgl(nijx_l), hfgb(nijx_l))
      ALLOCATE(hflr(nsolid,nijx_l), hflt(nsolid,nijx_l))
      ALLOCATE(hfll(nsolid,nijx_l), hflb(nsolid,nijx_l))
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
      IF (iturb .GE. 1) THEN
        kapgt = kapgt + kapg
      ELSE
        kapgt = kapg
      END IF
!
      CALL data_exchange(sieg)
      CALL data_exchange(siek)
      CALL data_exchange(tg)
      CALL data_exchange(tk)
      CALL data_exchange(kapgt)
!
      DO ij = 1, nij_l
        imesh = myij(0, 0, ij)
        IF(fl_l(ij).EQ.1) THEN
          CALL subscl(ij)
          CALL fsc_rt(egfr(ij), egft(ij), nb(sieg,ij),        &
               nb(rgp,ij), rnb(ug,ij), rnb(vg,ij), ij)
          CALL hotcg(hfgr(ij), hfgt(ij), hfgl(ij), hfgb(ij),  &
               nb(ep,ij), nb(tg,ij), nb(kapgt,ij), ij)
!
          egfr(ij) = egfr(ij) - hfgr(ij)
          egft(ij) = egft(ij) - hfgt(ij)
!
          DO k=1, nsolid
            CALL fsc_rt(elfr(k, ij), elft(k, ij), nb(siek(k,:),ij),     &
               nb(rlk(k,:),ij), rnb(uk(k,:),ij), rnb(vk(k,:),ij),ij)
            CALL hotck(hflr(k,ij), hflt(k,ij), hfll(k,ij), hflb(k,ij),  &
                 nb(rlk(k,:),ij), nb(tk(k,:),ij), ij, k)
!
          elfr(k, ij) = elfr(k, ij) - hflr(k, ij)
          elft(k, ij) = elft(k, ij) - hflt(k, ij)
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
      DO ij = 1, nij_l
        imesh = myij(0, 0, ij)
        IF(fl_l(ij).EQ.1) THEN
          CALL subscl(ij)
          j = ( imesh - 1 ) / nr + 1
          i = MOD( ( imesh - 1 ), nr) + 1
! 
          egfl(ij) = egfr(imj)
          egfb(ij) = egft(ijm)
          CALL fsc_lb(egfl(ij), egfb(ij), nb(sieg,ij),          &
               nb(rgp,ij), rnb(ug,ij), rnb(vg,ij), ij)
          CALL hotcg(hfgr(ij), hfgt(ij), hfgl(ij), hfgb(ij),    &
               nb(ep,ij), nb(tg,ij), nb(kapgt,ij), ij)
!
            egfl(ij) = egfl(ij) - hfgl(ij)
            egfb(ij) = egfb(ij) - hfgb(ij)
!
          DO k=1, nsolid
           IF (rlk(k,imj) * inrl(k) .LE. 1.D-9) THEN
            elfl(k, ij) = 0.0D0
           ELSE
            elfl(k, ij) = elfr(k, imj)
           END IF
           IF (rlk(k,ijm) * inrl(k) .LE. 1.D-9) THEN
            elfb(k, ij) = 0.0D0
           ELSE
            elfb(k, ij) = elft(k, ijm)
           END IF
            CALL fsc_lb(elfl(k, ij), elfb(k, ij), nb(siek(k,:),ij),      &
               nb(rlk(k,:),ij), rnb(uk(k,:),ij), rnb(vk(k,:),ij),ij)
            CALL hotck(hflr(k,ij), hflt(k,ij), hfll(k,ij), hflb(k,ij),   &
                 nb(rlk(k,:),ij), nb(tk(k,:),ij), ij, k)
!
          elfl(k, ij) = elfl(k, ij) - hfll(k, ij)
          elfb(k, ij) = elfb(k, ij) - hflb(k, ij)
!
          END DO
!
          c1 = dt * inr(i)*indr(i) * (egfr(ij)-egfl(ij))     &
             + dt * indz(j)        * (egft(ij)-egfb(ij))
!
          drc = dr(i)+(dr(i+1)+dr(i-1))/2.D0 
          dzc = dz(j)+(dz(j+1)+dz(j-1))/2.D0
          ugb = (ug(ij)+ug(imj))/2.D0
          vgb = (vg(ij)+vg(ijm))/2.D0
          upxy= dt/drc * ugb * (p(ijr)-p(ijl)) + dt/dzc * vgb * (p(ijt)-p(ijb))
!
          c2 = ep(ij) * (p(ij)-pn(ij)+upxy)
!
          rhg(ij) = - c1 + c2
!
          DO k=1,nsolid
            c3 = dt * inr(i)*indr(i) * (elfr(k,ij) - elfl(k,ij))    &
               + dt * indz(j)        * (elft(k,ij) - elfb(k,ij))
!
            rhk(k, ij) = - c3
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
