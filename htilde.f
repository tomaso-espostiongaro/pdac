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
! ... Compute the convective and diffusive enthalpy fluxes.
! ... (2D/3D-Compliant)
!
      USE control_flags, ONLY: job_type
      USE dimensions
      USE domain_decomposition, ONLY: ncint, ncdom, meshinds
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: sieg, sies, tg, ts
      USE gas_solid_temperature, ONLY: siegn, siesn
      USE grid, ONLY: dx, dy, dz, indx, indy, indz, inx
      USE grid, ONLY: fl_l
      USE particles_constants, ONLY: inrl
      USE pressure_epsilon, ONLY: ep
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm
      USE set_indexes, ONLY: ijke, ijkw, ijkn, ijks, ijkt, ijkb
      USE time_parameters, ONLY: dt
      IMPLICIT NONE
!
      REAL*8 :: egfx, egfy, egfz
      REAL*8 :: esfx, esfy, esfz
      REAL*8 :: hrexs, hrexg
      REAL*8 :: zero
      REAL*8 :: flx
      INTEGER :: is, m, l, is1
      INTEGER :: i, j, k, ijk, imesh
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
! ... Compute East, North, and Top fluxes in every cell
! ... within the computational domain
! 
      CALL compute_all_fluxes
!
! ... fluxes on left and bottom sides keep values
! ... entering from neighbouring cells.
!
      egfx = 0.D0; egfy = 0.D0; egfz = 0.D0
      esfx = 0.D0; esfy = 0.D0; esfz = 0.D0
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
          rhg(ijk) = - flx
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
      SUBROUTINE compute_all_fluxes

      USE control_flags, ONLY: job_type
      USE convective_fluxes_sc, ONLY: fsc, muscl_fsc
      USE diffusive_fluxes, ONLY: hotc
      USE dimensions, ONLY: nsolid
      USE domain_decomposition, ONLY: ncint, ncdom
      USE domain_decomposition, ONLY: data_exchange
      USE flux_limiters, ONLY: muscl
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: sieg, sies, tg, ts
      USE gas_solid_temperature, ONLY: siegn
      USE gas_solid_viscosity, ONLY: kapg
      USE gas_solid_viscosity, ONLY: gas_viscosity, part_viscosity
      USE grid, ONLY: fl_l
      USE particles_constants, ONLY: inrl, kap
      USE pressure_epsilon, ONLY: ep
      USE set_indexes, ONLY: stencil, cte
      USE set_indexes, ONLY: first_nb, first_rnb, third_nb, third_rnb
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes
      USE turbulence_model, ONLY: kapgt, iturb

      IMPLICIT NONE
!
      REAL*8 :: zero
      INTEGER :: is
      INTEGER :: ijk
      TYPE(stencil) :: u, v, w, dens, enth
      TYPE(stencil) :: eps, temp, kappa
!
! ... Exchange heat conductivity for diffusive fluxes
!
      CALL data_exchange(kapgt)

      IF (iturb >= 1) THEN
        kapgt = kapgt + kapg
      ELSE IF (iturb == 0) THEN
        kapgt = kapg
      END IF
!
      DO ijk = 1, ncint
        IF(fl_l(ijk) == 1) THEN
          CALL subscr(ijk)
!
! ... Here compute convective and diffusive fluxes (gas)
!
          IF (job_type == '2D') THEN
!
! ... First Order convective fluxes
!
            ! assemble first order computational stencils
            CALL first_nb(enth,sieg,ijk)
            CALL first_nb(dens,rgp,ijk)
            CALL first_rnb(u,ug,ijk)
            CALL first_rnb(w,wg,ijk)

            CALL fsc(egfe(ijk), egft(ijk),    &
                     egfe(imjk), egft(ijkm),    &
                     enth, dens, u, w, ijk)
!
! ... MUSCL second order correction
!
            IF (muscl > 0) THEN
              CALL third_nb(enth,sieg,ijk)
              CALL third_nb(dens,rgp,ijk)
              CALL muscl_fsc(egfe(ijk), egft(ijk),    & 
                             enth, dens, u, w, ijk)
            END IF
!
! ... Diffusive fluxes
!
            IF (gas_viscosity) THEN

              ! assemble first order computational stencils
              CALL first_nb(kappa,kapgt,ijk)
              CALL first_nb(eps,ep,ijk)
              CALL first_nb(temp,tg,ijk)
      
              CALL hotc(hgfe(ijk), zero, hgft(ijk),        &
                        hgfe(imjk), zero, hgft(ijkm),      &
                        eps, temp, kappa, ijk)
            END IF

          ELSE IF (job_type == '3D') THEN

! ... First Order convective fluxes
!
            ! assemble first order computational stencils
            CALL first_nb(enth,sieg,ijk)
            CALL first_nb(dens,rgp,ijk)
            CALL first_rnb(u,ug,ijk)
            CALL first_rnb(v,vg,ijk)
            CALL first_rnb(w,wg,ijk)

            CALL fsc(egfe(ijk), egfn(ijk), egft(ijk),    &
                     egfe(imjk), egfn(ijmk), egft(ijkm),    &
                     enth, dens, u, v, w, ijk)
!
! ... MUSCL second order correction
!
            IF (muscl > 0) THEN
              CALL third_nb(enth,sieg,ijk)
              CALL third_nb(dens,rgp,ijk)
              CALL muscl_fsc(egfe(ijk), egfn(ijk), egft(ijk),    &
                             enth, dens, u, v, w, ijk)
            END IF
!
! ... Diffusive fluxes
!
            IF (gas_viscosity) THEN

              ! assemble first order computational stencils
              CALL first_nb(kappa,kapgt,ijk)
              CALL first_nb(eps,ep,ijk)
              CALL first_nb(temp,tg,ijk)
      
              CALL hotc(hgfe(ijk), hgfn(ijk), hgft(ijk),         &
                        hgfe(imjk), hgfn(ijmk), hgft(ijkm),      &
                        eps, temp, kappa, ijk)
            END IF

          END IF
!
          IF (gas_viscosity) THEN
            egfe(ijk) = egfe(ijk) - hgfe(ijk)
            egft(ijk) = egft(ijk) - hgft(ijk)
!
            IF (fl_l(imjk) /= 1) egfe(imjk) = egfe(imjk) - hgfe(imjk)
            IF (fl_l(ijkm) /= 1) egft(ijkm) = egft(ijkm) - hgft(ijkm)

            IF (job_type == '3D') THEN
              egfn(ijk) = egfn(ijk) - hgfn(ijk)
              IF (fl_l(ijmk) /= 1) egfn(ijmk) = egfn(ijmk) - hgfn(ijmk)
            END IF
          END IF
!
! ... Here compute convective and diffusive fluxes (particles)
!
          DO is=1, nsolid
!
            IF (job_type == '2D') THEN

! ... First Order convective fluxes
!
              ! ... assemble first order computational stencils
              CALL first_nb(enth,sies(:,is),ijk)
              CALL first_nb(dens,rlk(:,is),ijk)
              CALL first_rnb(u, us(:,is),ijk)
              CALL first_rnb(w, ws(:,is),ijk)

              CALL fsc(esfe(ijk, is), esft(ijk, is),  &
                       esfe(imjk, is), esft(ijkm, is),  &
                       enth, dens, u, w, ijk)
!
! ... MUSCL second order correction
!
              IF (muscl > 0) THEN
                CALL third_nb(enth,sies(:,is),ijk)
                CALL third_nb(dens,rlk(:,is),ijk)
                CALL muscl_fsc(esfe(ijk, is), esft(ijk, is),  &
                               enth, dens, u, w, ijk)
              END IF
!
! ... Diffusive fluxes
!
              IF (part_viscosity) THEN

                ! ... assemble first order computational stencils
                eps   = inrl(is) * dens
                kappa = cte(kap(is))
                CALL first_nb(temp,ts(:,is),ijk)
            
                CALL hotc(hsfe(ijk,is), zero, hsft(ijk, is),    &
                          hsfe(imjk,is), zero, hsft(ijkm, is),    &
                          eps, temp, kappa, ijk)
              END IF

            ELSE IF (job_type == '3D') THEN

! ... First Order convective fluxes
!
              ! ... assemble first order computational stencils
              CALL first_nb(enth,sies(:,is),ijk)
              CALL first_nb(dens,rlk(:,is),ijk)
              CALL first_rnb(u, us(:,is),ijk)
              CALL first_rnb(v, vs(:,is),ijk)
              CALL first_rnb(w, ws(:,is),ijk)

              CALL fsc(esfe(ijk, is), esfn(ijk, is), esft(ijk, is),  &
                       esfe(imjk, is), esfn(ijmk, is), esft(ijkm, is),  &
                       enth, dens, u, v, w, ijk)
!
! ... MUSCL second order correction
!
              IF (muscl > 0) THEN
                CALL third_nb(enth,sies(:,is),ijk)
                CALL third_nb(dens,rlk(:,is),ijk)
                CALL muscl_fsc(esfe(ijk, is), esfn(ijk, is), esft(ijk, is),  &
                               enth, dens, u, v, w, ijk)
              END IF
!
! ... Diffusive fluxes
!
              IF (part_viscosity) THEN

                ! ... assemble first order computational stencils
                eps   = inrl(is) * dens
                kappa = cte(kap(is))
                CALL first_nb(temp,ts(:,is),ijk)
            
                CALL hotc(hsfe(ijk,is), hsfn(ijk, is), hsft(ijk, is),    &
                          hsfe(imjk,is), hsfn(ijmk, is), hsft(ijkm, is),    &
                          eps, temp, kappa, ijk)
              END IF

            END IF
!
            IF (part_viscosity) THEN
              esfe(ijk, is) = esfe(ijk, is) - hsfe(ijk,is)
              esft(ijk, is) = esft(ijk, is) - hsft(ijk, is)
!
              IF (fl_l(imjk) /= 1) esfe(imjk,is) = esfe(imjk,is) - hsfe(imjk,is)
              IF (fl_l(ijkm) /= 1) esft(ijkm,is) = esft(ijkm,is) - hsft(ijkm,is)

              IF (job_type == '3D') THEN
                esfn(ijk, is) = esfn(ijk, is) - hsfn(ijk, is)
                IF (fl_l(ijmk)/=1) esfn(ijmk,is)=esfn(ijmk,is)-hsfn(ijmk,is)
              END IF
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

      RETURN
      END SUBROUTINE compute_all_fluxes
!----------------------------------------------------------------------
      END MODULE tilde_energy
!-----------------------------------------------------------
