!-----------------------------------------------------------------------
      MODULE tilde_energy
!-----------------------------------------------------------------------
      IMPLICIT NONE
!
! ... tilde enthalpy terms
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: rhg
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rhs
!
! ... convective enthalpy fluxes
!
      REAL*8, ALLOCATABLE :: egfe(:), egfn(:), egft(:)
      REAL*8, ALLOCATABLE :: esfe(:,:), esfn(:,:), esft(:,:)
!
! ... diffusive enthalpy fluxes
!
      REAL*8, ALLOCATABLE :: hgfe(:), hgfn(:), hgft(:)
      REAL*8, ALLOCATABLE :: hsfe(:,:), hsfn(:,:), hsft(:,:)
!
      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE htilde
! ... (2D/3D-Compliant)
!
      USE convective_fluxes_sc, ONLY: fsc
      USE control_flags, ONLY: job_type
      USE diffusive_fluxes, ONLY: hotc
      USE dimensions
      USE domain_decomposition, ONLY: ncint, ncdom, myijk, data_exchange, meshinds
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: sieg, sies, tg, ts
      USE gas_solid_viscosity, ONLY: kapg
      USE grid, ONLY: dx, dy, dz, indx, indy, indz, inx
      USE grid, ONLY: fl_l
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE particles_constants, ONLY: inrl, kap
      USE pressure_epsilon, ONLY: p, pn, ep
      USE set_indexes
      USE time_parameters, ONLY: dt
      USE turbulence_model, ONLY: kapgt, iturb
      IMPLICIT NONE
!
      REAL*8 :: egfx, egfy, egfz
      REAL*8 :: esfx, esfy, esfz
      REAL*8 :: hrexs, hrexg
      REAL*8 :: zero
      REAL*8 :: flx, deltap, dpxyz
      REAL*8 :: indxc, indyc, indzc, ugc, vgc, wgc
      INTEGER :: is, m, l, is1
      INTEGER :: i, j, k, ijk, imesh
      TYPE(stencil) :: u, v, w, dens, enth
      TYPE(stencil) :: eps, temp, kappa
!
      ALLOCATE(rhg(ncint))
      ALLOCATE(rhs(ncint,nsolid))
!
      ALLOCATE(egfe(ncdom), egft(ncdom))
      ALLOCATE(esfe(ncdom,nsolid), esft(ncdom,nsolid))
!
      ALLOCATE(hgfe(ncdom), hgft(ncdom))
      ALLOCATE(hsfe(ncdom,nsolid), hsft(ncdom,nsolid))
!
      egfe = 0.0D0;  hgfe = 0.0D0
      egft = 0.0D0;  hgft = 0.0D0
      esfe = 0.0D0;  hsfe = 0.0D0
      esft = 0.0D0;  hsft = 0.0D0
!
      IF (job_type == '3D') THEN
        ALLOCATE(egfn(ncdom))
        ALLOCATE(esfn(ncdom,nsolid))
!
        ALLOCATE(hgfn(ncdom))
        ALLOCATE(hsfn(ncdom,nsolid))
!
        egfn = 0.0D0;  hgfn = 0.0D0
        esfn = 0.0D0;  hsfn = 0.0D0
      END IF
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
        IF(fl_l(ijk) == 1) THEN
          CALL subscr(ijk)
!
! ... compute convective and diffusive fluxes (gas)
!
! ... assemble the computational stencils
          CALL nb(enth,sieg,ijk)
          CALL nb(dens,rgp,ijk)
          CALL nb(eps,ep,ijk)
          CALL nb(temp,tg,ijk)
          CALL nb(kappa,kapgt,ijk)
          CALL rnb(u,ug,ijk)
          CALL rnb(w,wg,ijk)

          IF (job_type == '2D') THEN

            CALL fsc(egfe(ijk), egft(ijk),    &
                     egfe(imjk), egft(ijkm),    &
                     enth, dens, u, w, ijk)
!
            CALL hotc(hgfe(ijk), zero, hgft(ijk),      &
                      hgfe(imjk), zero, hgft(ijkm),      &
                      eps, temp, kappa, ijk)

          ELSE IF (job_type == '3D') THEN

            CALL rnb(v,vg,ijk)
            CALL fsc(egfe(ijk), egfn(ijk), egft(ijk),    &
                     egfe(imjk), egfn(ijmk), egft(ijkm),    &
                     enth, dens, u, v, w, ijk)
!
            CALL hotc(hgfe(ijk), hgfn(ijk), hgft(ijk),      &
                      hgfe(imjk), hgfn(ijmk), hgft(ijkm),      &
                      eps, temp, kappa, ijk)

          END IF
!
            egfe(ijk) = egfe(ijk) - hgfe(ijk)
            egft(ijk) = egft(ijk) - hgft(ijk)
!
            IF (fl_l(imjk) /= 1) egfe(imjk) = egfe(imjk) - hgfe(imjk)
            IF (fl_l(ijkm) /= 1) egft(ijkm) = egft(ijkm) - hgft(ijkm)

            IF (job_type == '3D') THEN
              egfn(ijk) = egfn(ijk) - hgfn(ijk)
              IF (fl_l(ijmk) /= 1) egfn(ijmk) = egfn(ijmk) - hgfn(ijmk)
            END IF
!
          DO is=1, nsolid
!
! ... convective and diffusive fluxes (particles)
!
! ... assemble computational stencils
            CALL nb(enth,sies(:,is),ijk)
            CALL nb(dens,rlk(:,is),ijk)
            CALL nb(temp,ts(:,is),ijk)
            CALL rnb(u, us(:,is),ijk)
            CALL rnb(w, ws(:,is),ijk)
            eps   = inrl(is) * dens
            kappa = cte(kap(is))
            
            IF (job_type == '2D') THEN

              CALL fsc(esfe(ijk, is), esft(ijk, is),  &
                       esfe(imjk, is), esft(ijkm, is),  &
                       enth, dens, u, w, ijk)
!
              CALL hotc(hsfe(ijk,is), zero, hsft(ijk, is),    &
                        hsfe(imjk,is), zero, hsft(ijkm, is),    &
                        eps, temp, kappa, ijk)

            ELSE IF (job_type == '3D') THEN

              CALL rnb(v, vs(:,is),ijk)
              CALL fsc(esfe(ijk, is), esfn(ijk, is), esft(ijk, is),  &
                       esfe(imjk, is), esfn(ijmk, is), esft(ijkm, is),  &
                       enth, dens, u, v, w, ijk)
!
              CALL hotc(hsfe(ijk,is), hsfn(ijk, is), hsft(ijk, is),    &
                        hsfe(imjk,is), hsfn(ijmk, is), hsft(ijkm, is),    &
                        eps, temp, kappa, ijk)

            END IF
!
            esfe(ijk, is) = esfe(ijk, is) - hsfe(ijk,is)
            esft(ijk, is) = esft(ijk, is) - hsft(ijk, is)
!
            IF (fl_l(imjk) /= 1) esfe(imjk,is) = esfe(imjk,is) - hsfe(imjk,is)
            IF (fl_l(ijkm) /= 1) esft(ijkm,is) = esft(ijkm,is) - hsft(ijkm,is)

            IF (job_type == '3D') THEN
              esfn(ijk, is) = esfn(ijk, is) - hsfn(ijk, is)
              IF (fl_l(ijmk) /= 1) esfn(ijmk,is) = esfn(ijmk,is) - hsfn(ijmk,is)
            END IF
!
          END DO
        END IF
      END DO
!
      CALL data_exchange(egfe)
      CALL data_exchange(egft)
      CALL data_exchange(esfe)
      CALL data_exchange(esft)

      IF (job_type == '3D') THEN
        CALL data_exchange(egfn)
        CALL data_exchange(esfn)
      END IF
!
! ... fluxes on left and bottom sides keep values
! ... entering from neighbouring cells.
!
      egfx = 0.D0; egfy = 0.D0; egfz = 0.D0
      esfx = 0.D0; esfy = 0.D0; esfz = 0.D0
      indxc = 0.D0 ; indyc = 0.D0 ; indzc = 0.D0
      ugc = 0.D0 ; vgc = 0.D0 ; wgc = 0.D0
!
      DO ijk = 1, ncint
        IF(fl_l(ijk) == 1) THEN
          CALL meshinds(ijk,imesh,i,j,k)
          CALL subscr(ijk)
! 
          egfx = egfe(ijk) - egfe(imjk)
          egfz = egft(ijk) - egft(ijkm)
!
          IF (job_type == '3D') THEN
            egfy = egfn(ijk) - egfn(ijmk)
          END IF

          flx = dt * indx(i) * egfx * inx(i) +   &
                dt * indy(j) * egfy          +   &
                dt * indz(k) * egfz
!
          indxc = 1.D0/(dx(i)+(dx(i+1)+dx(i-1))*0.5D0)
          indzc = 1.D0/(dz(k)+(dz(k+1)+dz(k-1))*0.5D0)
          ugc = (ug(ijk)+ug(imjk))/2.D0
          wgc = (wg(ijk)+wg(ijkm))/2.D0

          IF (job_type == '3D') THEN
            indyc = 1.D0 / (dy(j)+(dy(j+1)+dy(j-1))*0.5D0)
            vgc = (vg(ijk)+vg(ijmk))/2.D0
          END IF
!
          dpxyz= dt * indxc * ugc * (p(ijke)-p(ijkw)) +   &
                 dt * indyc * vgc * (p(ijkn)-p(ijks)) +   &
                 dt * indzc * wgc * (p(ijkt)-p(ijkb))
!
          deltap = ep(ijk) * (p(ijk) - pn(ijk) + dpxyz)
!
          rhg(ijk) = - flx + deltap
!
! ... Same procedure carried out for solids
!
          DO is=1, nsolid
           IF (rlk(imjk,is) * inrl(is) <= 1.D-9) THEN
             esfx = esfe(ijk, is)
           ELSE
             esfx = esfe(ijk, is) - esfe(imjk, is)
           END IF

           IF (rlk(ijkm,is) * inrl(is) <= 1.D-9) THEN
             esfz = esft(ijk, is)
           ELSE
             esfz = esft(ijk, is) - esft(ijkm, is)
           END IF

           IF (job_type == '3D') THEN
             IF (rlk(ijmk,is) * inrl(is) <= 1.D-9) THEN
               esfy = esfn(ijk, is)
             ELSE
               esfy = esfn(ijk, is) - esfn(ijmk, is)
             END IF
           END IF
!
            flx = dt * indx(i) * esfx * inx(i) +  &
                  dt * indy(j) * esfy          +  &
                  dt * indz(k) * esfz
!
            rhs(ijk,is) = - flx
          END DO
!
        END IF
      END DO
!
!      CALL test_fluxes
!
      DEALLOCATE(egfe, egft)
      DEALLOCATE(esfe, esft)
!
      DEALLOCATE(hgfe, hgft)
      DEALLOCATE(hsfe, hsft)
!
      IF (job_type == '3D') THEN
        DEALLOCATE(egfn)
        DEALLOCATE(esfn)
        DEALLOCATE(hgfn)
        DEALLOCATE(hsfn)
      END IF
!
      RETURN
      END SUBROUTINE htilde
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
      filnam='output.htest'

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
        CALL write_array( 12, rhg, sgl, lform )
        CALL write_array( 12, egfe, sgl, lform )
        CALL write_array( 12, egft, sgl, lform )
        CALL write_array( 12, hgfe, sgl, lform )
        CALL write_array( 12, hgft, sgl, lform )
      ELSE IF (job_type == '3D') THEN
        CALL write_array( 12, rhg, sgl, lform )
        CALL write_array( 12, egfe, sgl, lform )
        CALL write_array( 12, egfn, sgl, lform )
        CALL write_array( 12, egft, sgl, lform )
        CALL write_array( 12, hgfe, sgl, lform )
        CALL write_array( 12, hgfn, sgl, lform )
        CALL write_array( 12, hgft, sgl, lform )
      END IF
!
      DO is = 1, nsolid
        IF (job_type == '2D') THEN
          CALL write_array( 12, rhs(:,is), sgl, lform )
          CALL write_array( 12, esfe(:,is), sgl, lform )
          CALL write_array( 12, esft(:,is), sgl, lform )
          CALL write_array( 12, hsfe(:,is), sgl, lform )
          CALL write_array( 12, hsft(:,is), sgl, lform )
        ELSE IF (job_type == '3D') THEN
          CALL write_array( 12, rhs(:,is), sgl, lform )
          CALL write_array( 12, esfe(:,is), sgl, lform )
          CALL write_array( 12, esfn(:,is), sgl, lform )
          CALL write_array( 12, esft(:,is), sgl, lform )
          CALL write_array( 12, hsfe(:,is), sgl, lform )
          CALL write_array( 12, hsfn(:,is), sgl, lform )
          CALL write_array( 12, hsft(:,is), sgl, lform )
        END IF
      END DO

      IF( mpime == root ) THEN
        CLOSE (12)
      END IF
!
      RETURN
      END SUBROUTINE test_fluxes
!-----------------------------------------------------------
      END MODULE tilde_energy
!-----------------------------------------------------------
