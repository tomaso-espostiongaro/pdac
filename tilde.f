!----------------------------------------------------------------------
      MODULE tilde_momentum
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: rug, rwg
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rus, rws
      REAL*8, DIMENSION(:),   ALLOCATABLE :: rugn, rwgn
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rusn, rwsn
!
      REAL*8, ALLOCATABLE :: ugfr(:), ugft(:)
      REAL*8, ALLOCATABLE :: ugfl(:), ugfb(:)
      REAL*8, ALLOCATABLE :: wgfr(:), wgft(:)
      REAL*8, ALLOCATABLE :: wgfl(:), wgfb(:)
!
      REAL*8, ALLOCATABLE :: ulfr(:,:), ulft(:,:)
      REAL*8, ALLOCATABLE :: ulfl(:,:), ulfb(:,:)
      REAL*8, ALLOCATABLE :: wlfr(:,:), wlft(:,:)
      REAL*8, ALLOCATABLE :: wlfl(:,:), wlfb(:,:)
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE  :: kpgv
      REAL*8, DIMENSION(:,:), ALLOCATABLE  :: appu, appw
!
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------

! ... MODIFICARE_X3D


!----------------------------------------------------------------------
      SUBROUTINE euvel
!----------------------------------------------------------------------
! .. This routine computes explicitly gas and particles momenta at time ndt
!
      USE dimensions
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_velocity, ONLY: ug, wg, us, ws
      USE grid, ONLY: dz, dr, fl_l
      USE grid, ONLY: inr, inrb, indr, indz
      USE grid, ONLY: nij_l, myij, nijx_l, data_exchange
      USE set_indexes
!
      IMPLICIT NONE
      SAVE
!
      REAL*8 :: drp, dzp, indrp, indzp
      REAL*8 :: rgp_e, rgp_n, rlk_e, rlk_n
      INTEGER :: is, i, j, ij, imesh
!
      IF (ALLOCATED(rugn)) DEALLOCATE(rugn)
      IF (ALLOCATED(rwgn)) DEALLOCATE(rwgn)
      IF (ALLOCATED(rusn)) DEALLOCATE(rusn)
      IF (ALLOCATED(rwsn)) DEALLOCATE(rwsn)
      ALLOCATE( rugn(nij_l),  rwgn(nij_l))
      ALLOCATE( rusn(nsolid,nij_l),  rwsn(nsolid,nij_l))
      rugn = 0.0
      rwgn = 0.0
      rusn = 0.0
      rwsn = 0.0
!
      CALL data_exchange(rgp)
      CALL data_exchange(rlk)
!
      DO ij = 1, nij_l
       imesh = myij(0,0,ij)
       IF(fl_l(ij).EQ.1) THEN
         CALL subscr(ij)
         j = ( imesh - 1 ) / nr + 1
         i = MOD( ( imesh - 1 ), nr) + 1
!
         drp=dr(i)+dr(i+1)
         dzp=dz(j)+dz(j+1)
         indrp=1.D0/drp
         indzp=1.D0/dzp
!
         rgp_e = (dr(i+1)*rgp(ij)+dr(i)*rgp(ijr))*indrp
         rgp_n = (dz(j+1)*rgp(ij)+dz(j)*rgp(ijt))*indzp
         rugn(ij) = rgp_e * ug(ij)
         rwgn(ij) = rgp_n * wg(ij)
!
         DO is=1, nsolid
          rlk_e = (rlk(is,ij)*dr(i+1)+rlk(is,ijr)*dr(i))*indrp
          rlk_n = (rlk(is,ij)*dz(j+1)+rlk(is,ijt)*dz(j))*indzp
          rusn(is,ij) = rlk_e * us(is,ij)
          rwsn(is,ij) = rlk_n * ws(is,ij)
         END DO
        END IF
      END DO
!
      RETURN
      END SUBROUTINE euvel
!----------------------------------------------------------------------
      SUBROUTINE tilde
!
      USE atmosphere, ONLY: gravx, gravz
      USE grid, ONLY: fl_l
      USE dimensions
      USE eulerian_flux, ONLY: fu_lb, fu_rt, fw_lb, fw_rt
      USE gas_solid_density, ONLY: rog, rgp, rlk
      USE gas_solid_velocity, ONLY: ug, wg, us, ws
      USE grid, ONLY: dz, dr
      USE grid, ONLY: inr, inrb, indr, indz
      USE momentum_transfer, ONLY: kdrags, inter
      USE particles_constants, ONLY: phi, epsl, dkf, epsu, dk, rl, inrl, philim
      USE pressure_epsilon, ONLY: ep
      USE time_parameters, ONLY: dt
      USE turbulence, ONLY: iss, iturb
      USE gas_solid_viscosity, ONLY: viscg, viscs
      USE gas_solid_viscosity, ONLY: mug
      USE gas_solid_viscosity, ONLY: gvisx, gvisz, pvisx, pvisz
      USE grid, ONLY: nij_l, myij, nijx_l, data_exchange
      USE parallel, ONLY: mpime
      USE set_indexes
!
      IMPLICIT NONE
      SAVE
!
      INTEGER :: i, j, is, imesh
      INTEGER :: ij
      REAL*8 :: drp, dzp, indrp, indzp
      REAL*8 :: ugfx, ugfz, wgfx, wgfz
      REAL*8 :: ulfx, ulfz, wlfx, wlfz
!
      ALLOCATE(gvisx(nij_l), gvisz(nij_l))
      ALLOCATE(pvisx(nsolid, nij_l), pvisz(nsolid, nij_l))
      gvisx = 0.D0; gvisz = 0.D0
      pvisx = 0.D0; pvisz = 0.D0

      IF (iturb .LT. 1) THEN
        CALL data_exchange(ug)
        CALL data_exchange(wg)
      END IF
      IF (iss .LT. 1) THEN
        CALL data_exchange(us)
        CALL data_exchange(ws)
      END IF
      CALL data_exchange(ep)
      CALL data_exchange(rgp)
      CALL data_exchange(rlk)
!
! ... Calculate gvisx gvisz (gas viscous stress tensor).
!
      CALL viscg        
!
! ... Calculate pvisx pvisz (particles viscous stress tensor).
!
      DO is=1,nsolid
        CALL viscs(is)   
      END DO
!
! ... Allocate and initialize local arrays (gas).
!
      ALLOCATE( rug(nijx_l),  rwg(nijx_l))
      ALLOCATE(ugfr(nijx_l), ugft(nijx_l))
      ALLOCATE(ugfl(nijx_l), ugfb(nijx_l))
      ALLOCATE(wgfr(nijx_l), wgft(nijx_l))
      ALLOCATE(wgfl(nijx_l), wgfb(nijx_l))

      rug = 0.0D0
      rwg = 0.0D0
      ugfr = 0.0D0
      ugft = 0.0D0
      ugfl = 0.0D0
      ugfb = 0.0D0
      wgfr = 0.0D0
      wgft = 0.0D0
      wgfl = 0.0D0
      wgfb = 0.0D0
!
! ... Allocate and initialize local arrays (particles).
!
      ALLOCATE(rus(nsolid, nijx_l), rws(nsolid, nijx_l) )
      ALLOCATE(ulfr(nsolid, nijx_l), ulft(nsolid, nijx_l))
      ALLOCATE(ulfl(nsolid, nijx_l), ulfb(nsolid, nijx_l))
      ALLOCATE(wlfr(nsolid, nijx_l), wlft(nsolid, nijx_l))
      ALLOCATE(wlfl(nsolid, nijx_l), wlfb(nsolid, nijx_l))

      rus  = 0.0D0
      rws  = 0.0D0
      ulfr = 0.0D0
      ulft = 0.0D0
      ulfl = 0.0D0
      ulfb = 0.0D0
      wlfr = 0.0D0
      wlft = 0.0D0
      wlfl = 0.0D0
      wlfb = 0.0D0
!
! ... Allocate and initialize interphase terms.
!
      ALLOCATE(kpgv(nsolid,nij_l))
      ALLOCATE(appu(((nsolid+1)**2+(nsolid+1))/2, nijx_l),   &
               appw(((nsolid+1)**2+(nsolid+1))/2, nijx_l))

      kpgv = 0.0D0
      appu = 0.0D0
      appw = 0.0D0
! 
! ... Compute fluxes on right and top sides of a cell
! ... in the whole computational domain.
!
      DO ij = 1, nij_l
        imesh = myij(0,0,ij)
        IF(fl_l(ij).EQ.1) THEN
          CALL subscr(ij)
          CALL fu_rt(ugfr(ij), ugft(ij), nb(rgp,ij),             &
                     rnb(ug,ij), rnb(wg,ij), ij)
          CALL fw_rt(wgfr(ij), wgft(ij), nb(rgp,ij),             &
                     rnb(ug,ij), rnb(wg,ij), ij)
!
          DO is = 1, nsolid
            CALL fu_rt(ulfr(is,ij), ulft(is,ij), nb(rlk(is,:),ij),  &
                       rnb(us(is,:),ij), rnb(ws(is,:),ij), ij)
            CALL fw_rt(wlfr(is,ij), wlft(is,ij), nb(rlk(is,:),ij),  &
                       rnb(us(is,:),ij), rnb(ws(is,:),ij), ij)
          END DO
        END IF         
      END DO
!
      CALL data_exchange(ugfr)
      CALL data_exchange(ugft)
      CALL data_exchange(wgfr)
      CALL data_exchange(wgft)
      CALL data_exchange(ulfr)
      CALL data_exchange(ulft)
      CALL data_exchange(wlfr)
      CALL data_exchange(wlft)
!
! ... fluxes on left and bottom sides keep values 
! ... of right and top fluxes from neighbouring (left and bottom) cells.
! ... On boundaries, fluxes on left and bottom sides must be calculated
!
      DO ij = 1, nij_l
        imesh = myij(0,0,ij)
        IF(fl_l(ij).EQ.1) THEN
          CALL subscr(ij)
          j = ( imesh - 1 ) / nr + 1
          i = MOD( ( imesh - 1 ), nr) + 1
!
          drp=dr(i)+dr(i+1)
          dzp=dz(j)+dz(j+1)
          indrp=1.D0/drp
          indzp=1.D0/dzp
!
! ... left and bottom fluxes (gas)
!
          ugfl(ij) = ugfr(imj)
          ugfb(ij) = ugft(ijm)
          wgfl(ij) = wgfr(imj) 
          wgfb(ij) = wgft(ijm)
!
          CALL fu_lb(ugfl(ij), ugfb(ij), nb(rgp,ij), rnb(ug,ij), rnb(wg,ij), ij)
          CALL fw_lb(wgfl(ij), wgfb(ij), nb(rgp,ij), rnb(ug,ij), rnb(wg,ij), ij)
!
! ... compute the flux balance in the radial (x)
! ... and vertical (z) directions, for the gas phase
!
          ugfx = ugfr(ij) - ugfl(ij)
          ugfz = ugft(ij) - ugfb(ij)
          wgfx = wgfr(ij) - wgfl(ij)
          wgfz = wgft(ij) - wgfb(ij)
!
! ... compute explicit (tilde) terms in the momentum equation (gas)
! 
          rug(ij) = rugn(ij)                                         &
     &      + dt * (dr(i+1)*rgp(ij)+dr(i)*rgp(ijr))*indrp * gravx    &
     &      - dt * inrb(i) * indrp * 2.D0 * ugfx                     &
     &      - dt * indz(j) * ugfz                                    &
     &      + dt * gvisx(ij)                  

          rwg(ij) = rwgn(ij)                                         &
     &      + dt * (dz(j+1)*rgp(ij)+dz(j)*rgp(ijt))*indzp * gravz    &
     &      - dt * inr(i) * indr(i) * wgfx                           &
     &      - dt * indzp * 2.D0 * wgfz                               &
     &      + dt * gvisz(ij)
!
! ... same procedure carried out for particulate phases
!
          DO is = 1, nsolid
!
! ... left and bottom fluxes (particles)
! 
            ulfl(is,ij) = ulfr(is,imj)
            ulfb(is,ij) = ulft(is,ijm)
            wlfl(is,ij) = wlfr(is,imj)
            wlfb(is,ij) = wlft(is,ijm)
!
            CALL fu_lb(ulfl(is,ij), ulfb(is,ij), nb(rlk(is,:),ij),      &
     &           rnb(us(is,:),ij), rnb(ws(is,:),ij), ij)
            CALL fw_lb(wlfl(is,ij), wlfb(is,ij), nb(rlk(is,:),ij),      &
     &           rnb(us(is,:),ij), rnb(ws(is,:),ij), ij)
!
! ... compute explicit (tilde) terms in the momentum equation (particles)
! 
              ulfx = ulfr(is,ij) - ulfl(is,ij)
              ulfz = ulft(is,ij) - ulfb(is,ij)

              rus(is,ij) = rusn(is,ij)                                 &
     &         + dt*(rlk(is,ij)*dr(i+1)+rlk(is,ijr)*dr(i))*indrp*gravx &
     &         - dt*inrb(i)*indrp*2.D0*ulfx                          &
     &         - dt*indz(j)*ulfz                                     &
     &         + dt*pvisx(is,ij)

              wlfx = wlfr(is,ij) - wlfl(is,ij)
              wlfz = wlft(is,ij) - wlfb(is,ij)

              rws(is,ij) = rwsn(is,ij)                                 &
     &         + dt*(rlk(is,ij)*dz(j+1)+rlk(is,ijt)*dz(j))*indzp*gravz &
     &         - dt*inr(i)*indr(i)* wlfx                             &
     &         - dt*indzp*2.D0* wlfz                                 &
     &         + dt*pvisz(is,ij)    

          END DO 
!
! ... Compute interphase coefficients
!
          CALL kdrags(kpgv(:,ij), ug(ij), ug(imj),                          &
     &                wg(ij), wg(ijm), us(:,ij), us(:,imj),                 &
     &                ws(:,ij), ws(:,ijm), ep(ij),                          &
     &                rog(ij), rgp(ij), rlk(:,ij), mug(ij))                  
          CALL inter(appu(:,ij), appw(:,ij), kpgv(:,ij),                    &
     &               us(:,ij), us(:,imj), ws(:,ij), ws(:,ijm), rlk(:,ij))
!
        END IF
      END DO
!
      DEALLOCATE(ugfr, ugft)
      DEALLOCATE(ugfl, ugfb)
      DEALLOCATE(wgfr, wgft)
      DEALLOCATE(wgfl, wgfb)
      DEALLOCATE(gvisx, gvisz)
!
      DEALLOCATE(ulfr, ulft)
      DEALLOCATE(ulfl, ulfb)
      DEALLOCATE(wlfr, wlft)
      DEALLOCATE(wlfl, wlfb)
      DEALLOCATE(pvisx, pvisz)
!
      DEALLOCATE(kpgv)
!
      CALL data_exchange(appu)
      CALL data_exchange(appw)
      CALL data_exchange(rug)
      CALL data_exchange(rwg)
      CALL data_exchange(rus)
      CALL data_exchange(rws)
!
      RETURN
      END SUBROUTINE tilde
!----------------------------------------------------------------------
      END MODULE tilde_momentum
!----------------------------------------------------------------------
