!----------------------------------------------------------------------
      MODULE tilde_momentum
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: rug, rvg
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: ruk, rvk
      REAL*8, DIMENSION(:),   ALLOCATABLE :: rugn, rvgn
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rukn, rvkn
!
      REAL*8, ALLOCATABLE :: ugfr(:), ugft(:)
      REAL*8, ALLOCATABLE :: ugfl(:), ugfb(:)
      REAL*8, ALLOCATABLE :: vgfr(:), vgft(:)
      REAL*8, ALLOCATABLE :: vgfl(:), vgfb(:)
!
      REAL*8, ALLOCATABLE :: ulfr(:,:), ulft(:,:)
      REAL*8, ALLOCATABLE :: ulfl(:,:), ulfb(:,:)
      REAL*8, ALLOCATABLE :: vlfr(:,:), vlft(:,:)
      REAL*8, ALLOCATABLE :: vlfl(:,:), vlfb(:,:)
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE  :: kpgv
      REAL*8, DIMENSION(:,:), ALLOCATABLE  :: appu, appv
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
      USE gas_solid_velocity, ONLY: ug, vg, uk, vk
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
      INTEGER :: k, i, j, ij, imesh
!
      IF (ALLOCATED(rugn)) DEALLOCATE(rugn)
      IF (ALLOCATED(rvgn)) DEALLOCATE(rvgn)
      IF (ALLOCATED(rukn)) DEALLOCATE(rukn)
      IF (ALLOCATED(rvkn)) DEALLOCATE(rvkn)
      ALLOCATE( rugn(nij_l),  rvgn(nij_l))
      ALLOCATE( rukn(nsolid,nij_l),  rvkn(nsolid,nij_l))
      rugn = 0.0
      rvgn = 0.0
      rukn = 0.0
      rvkn = 0.0
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
         rvgn(ij) = rgp_n * vg(ij)
!
         DO k=1, nsolid
          rlk_e = (rlk(k,ij)*dr(i+1)+rlk(k,ijr)*dr(i))*indrp
          rlk_n = (rlk(k,ij)*dz(j+1)+rlk(k,ijt)*dz(j))*indzp
          rukn(k,ij) = rlk_e * uk(k,ij)
          rvkn(k,ij) = rlk_n * vk(k,ij)
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
      USE eulerian_flux, ONLY: fu_lb, fu_rt, fv_lb, fv_rt
      USE gas_solid_density, ONLY: rog, rgp, rlk
      USE gas_solid_velocity, ONLY: ug, vg, uk, vk
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
      INTEGER :: i, j, k, imesh
      INTEGER :: ij
      REAL*8 :: drp, dzp, indrp, indzp
      REAL*8 :: ugfx, ugfy, vgfx, vgfy
      REAL*8 :: ulfx, ulfy, vlfx, vlfy
!
      ALLOCATE(gvisx(nij_l), gvisz(nij_l))
      ALLOCATE(pvisx(nsolid, nij_l), pvisz(nsolid, nij_l))
      gvisx = 0.D0; gvisz = 0.D0
      pvisx = 0.D0; pvisz = 0.D0

      IF (iturb .LT. 1) THEN
        CALL data_exchange(ug)
        CALL data_exchange(vg)
      END IF
      IF (iss .LT. 1) THEN
        CALL data_exchange(uk)
        CALL data_exchange(vk)
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
      DO k=1,nsolid
        CALL viscs(k)   
      END DO
!
! ... Allocate and initialize local arrays (gas).
!
      ALLOCATE( rug(nijx_l),  rvg(nijx_l))
      ALLOCATE(ugfr(nijx_l), ugft(nijx_l))
      ALLOCATE(ugfl(nijx_l), ugfb(nijx_l))
      ALLOCATE(vgfr(nijx_l), vgft(nijx_l))
      ALLOCATE(vgfl(nijx_l), vgfb(nijx_l))

      rug = 0.0D0
      rvg = 0.0D0
      ugfr = 0.0D0
      ugft = 0.0D0
      ugfl = 0.0D0
      ugfb = 0.0D0
      vgfr = 0.0D0
      vgft = 0.0D0
      vgfl = 0.0D0
      vgfb = 0.0D0
!
! ... Allocate and initialize local arrays (particles).
!
      ALLOCATE(ruk(nsolid, nijx_l), rvk(nsolid, nijx_l) )
      ALLOCATE(ulfr(nsolid, nijx_l), ulft(nsolid, nijx_l))
      ALLOCATE(ulfl(nsolid, nijx_l), ulfb(nsolid, nijx_l))
      ALLOCATE(vlfr(nsolid, nijx_l), vlft(nsolid, nijx_l))
      ALLOCATE(vlfl(nsolid, nijx_l), vlfb(nsolid, nijx_l))

      ruk  = 0.0D0
      rvk  = 0.0D0
      ulfr = 0.0D0
      ulft = 0.0D0
      ulfl = 0.0D0
      ulfb = 0.0D0
      vlfr = 0.0D0
      vlft = 0.0D0
      vlfl = 0.0D0
      vlfb = 0.0D0
!
! ... Allocate and initialize interphase terms.
!
      ALLOCATE(kpgv(nsolid,nij_l))
      ALLOCATE(appu(((nsolid+1)**2+(nsolid+1))/2, nijx_l),   &
               appv(((nsolid+1)**2+(nsolid+1))/2, nijx_l))

      kpgv = 0.0D0
      appu = 0.0D0
      appv = 0.0D0
! 
! ... Compute fluxes on right and top sides of a cell
! ... in the whole computational domain.
!
      DO ij = 1, nij_l
        imesh = myij(0,0,ij)
        IF(fl_l(ij).EQ.1) THEN
          CALL subscr(ij)
          CALL fu_rt(ugfr(ij), ugft(ij), nb(rgp,ij),             &
                     rnb(ug,ij), rnb(vg,ij), ij)
          CALL fv_rt(vgfr(ij), vgft(ij), nb(rgp,ij),             &
                     rnb(ug,ij), rnb(vg,ij), ij)
!
          DO k = 1, nsolid
            CALL fu_rt(ulfr(k,ij), ulft(k,ij), nb(rlk(k,:),ij),  &
                       rnb(uk(k,:),ij), rnb(vk(k,:),ij), ij)
            CALL fv_rt(vlfr(k,ij), vlft(k,ij), nb(rlk(k,:),ij),  &
                       rnb(uk(k,:),ij), rnb(vk(k,:),ij), ij)
          END DO
        END IF         
      END DO
!
      CALL data_exchange(ugfr)
      CALL data_exchange(ugft)
      CALL data_exchange(vgfr)
      CALL data_exchange(vgft)
      CALL data_exchange(ulfr)
      CALL data_exchange(ulft)
      CALL data_exchange(vlfr)
      CALL data_exchange(vlft)
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
          vgfl(ij) = vgfr(imj) 
          vgfb(ij) = vgft(ijm)
!
          CALL fu_lb(ugfl(ij), ugfb(ij), nb(rgp,ij), rnb(ug,ij), rnb(vg,ij), ij)
          CALL fv_lb(vgfl(ij), vgfb(ij), nb(rgp,ij), rnb(ug,ij), rnb(vg,ij), ij)
!
! ... compute the flux balance in the radial (x)
! ... and vertical (y) directions, for the gas phase
!
          ugfx = ugfr(ij) - ugfl(ij)
          ugfy = ugft(ij) - ugfb(ij)
          vgfx = vgfr(ij) - vgfl(ij)
          vgfy = vgft(ij) - vgfb(ij)
!
! ... compute explicit (tilde) terms in the momentum equation (gas)
! 
          rug(ij) = rugn(ij)                                         &
     &      + dt * (dr(i+1)*rgp(ij)+dr(i)*rgp(ijr))*indrp * gravx    &
     &      - dt * inrb(i) * indrp * 2.D0 * ugfx                     &
     &      - dt * indz(j) * ugfy                                    &
     &      + dt * gvisx(ij)                  

          rvg(ij) = rvgn(ij)                                         &
     &      + dt * (dz(j+1)*rgp(ij)+dz(j)*rgp(ijt))*indzp * gravz    &
     &      - dt * inr(i) * indr(i) * vgfx                           &
     &      - dt * indzp * 2.D0 * vgfy                               &
     &      + dt * gvisz(ij)
!
! ... same procedure carried out for particulate phases
!
          DO k = 1, nsolid
!
! ... left and bottom fluxes (particles)
! 
            ulfl(k,ij) = ulfr(k,imj)
            ulfb(k,ij) = ulft(k,ijm)
            vlfl(k,ij) = vlfr(k,imj)
            vlfb(k,ij) = vlft(k,ijm)
!
            CALL fu_lb(ulfl(k,ij), ulfb(k,ij), nb(rlk(k,:),ij),      &
     &           rnb(uk(k,:),ij), rnb(vk(k,:),ij), ij)
            CALL fv_lb(vlfl(k,ij), vlfb(k,ij), nb(rlk(k,:),ij),      &
     &           rnb(uk(k,:),ij), rnb(vk(k,:),ij), ij)
!
! ... compute explicit (tilde) terms in the momentum equation (particles)
! 
              ulfx = ulfr(k,ij) - ulfl(k,ij)
              ulfy = ulft(k,ij) - ulfb(k,ij)

              ruk(k,ij) = rukn(k,ij)                                 &
     &         + dt*(rlk(k,ij)*dr(i+1)+rlk(k,ijr)*dr(i))*indrp*gravx &
     &         - dt*inrb(i)*indrp*2.D0*ulfx                          &
     &         - dt*indz(j)*ulfy                                     &
     &         + dt*pvisx(k,ij)

              vlfx = vlfr(k,ij) - vlfl(k,ij)
              vlfy = vlft(k,ij) - vlfb(k,ij)

              rvk(k,ij) = rvkn(k,ij)                                 &
     &         + dt*(rlk(k,ij)*dz(j+1)+rlk(k,ijt)*dz(j))*indzp*gravz &
     &         - dt*inr(i)*indr(i)* vlfx                             &
     &         - dt*indzp*2.D0* vlfy                                 &
     &         + dt*pvisz(k,ij)    

          END DO 
!
! ... Compute interphase coefficients
!
          CALL kdrags(kpgv(:,ij), ug(ij), ug(imj),                          &
     &                vg(ij), vg(ijm), uk(:,ij), uk(:,imj),                 &
     &                vk(:,ij), vk(:,ijm), ep(ij),                          &
     &                rog(ij), rgp(ij), rlk(:,ij), mug(ij))                  
          CALL inter(appu(:,ij), appv(:,ij), kpgv(:,ij),                    &
     &               uk(:,ij), uk(:,imj), vk(:,ij), vk(:,ijm), rlk(:,ij))
!
        END IF
      END DO
!
      DEALLOCATE(ugfr, ugft)
      DEALLOCATE(ugfl, ugfb)
      DEALLOCATE(vgfr, vgft)
      DEALLOCATE(vgfl, vgfb)
      DEALLOCATE(gvisx, gvisz)
!
      DEALLOCATE(ulfr, ulft)
      DEALLOCATE(ulfl, ulfb)
      DEALLOCATE(vlfr, vlft)
      DEALLOCATE(vlfl, vlfb)
      DEALLOCATE(pvisx, pvisz)
!
      DEALLOCATE(kpgv)
!
      CALL data_exchange(appu)
      CALL data_exchange(appv)
      CALL data_exchange(rug)
      CALL data_exchange(rvg)
      CALL data_exchange(ruk)
      CALL data_exchange(rvk)
!
      RETURN
      END SUBROUTINE tilde
!----------------------------------------------------------------------
      END MODULE tilde_momentum
!----------------------------------------------------------------------
