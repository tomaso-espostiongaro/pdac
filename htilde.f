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
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions
      USE domain_mapping, ONLY: ncint, ncdom, meshinds
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: sieg, sies, tg, ts
      USE gas_solid_temperature, ONLY: siegn, siesn
      USE grid, ONLY: dx, dy, dz, indx, indy, indz, inr
      USE grid, ONLY: flag
      USE immersed_boundaries, ONLY: faces, immb
      USE particles_constants, ONLY: inrl
      USE pressure_epsilon, ONLY: ep
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm
      USE set_indexes, ONLY: ctu1_subscr, ctu2_subscr, ctu3_subscr
      USE mass_sink, ONLY: sink
      USE time_parameters, ONLY: dt
      USE flux_limiters, ONLY: ctu
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
      INTEGER :: b_e, b_w, b_t, b_b, b_n, b_s
      REAL*8 :: ivf
      LOGICAL :: compute
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
      IF (job_type == JOB_TYPE_3D) THEN
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
      CALL compute_all_fluxes_htilde
!      CALL test_fluxes_htilde
!
! ... Fluxes on west, south, and bottom sides keep values
! ... entering from neighbouring cells.
!
      egfx = 0.D0; egfy = 0.D0; egfz = 0.D0
      esfx = 0.D0; esfy = 0.D0; esfz = 0.D0
!
      DO ijk = 1, ncint
        compute = BTEST(flag(ijk),0)
        IF( compute ) THEN

          ! ... Compute the volumes partially filled by the
          ! ... topography
          CALL faces(ijk, b_e, b_w, b_t, b_b, b_n, b_s, ivf)

          CALL meshinds(ijk,imesh,i,j,k)
          CALL subscr(ijk)
          IF (ctu > 0) CALL ctu1_subscr(ijk)
          IF (ctu > 1) CALL ctu2_subscr(ijk)
          IF (ctu > 2) CALL ctu3_subscr(ijk)
! 
          egfx = b_e * egfe(ijk) - b_w * egfe(imjk)
          egfz = b_t * egft(ijk) - b_b * egft(ijkm)
!
          IF (job_type == JOB_TYPE_3D) THEN
            egfy = b_n * egfn(ijk) - b_s * egfn(ijmk)
          END IF

          flx = dt * indx(i) * egfx * inr(i) +   &
                dt * indy(j) * egfy          +   &
                dt * indz(k) * egfz
          !
          ! ... volume correction for cells partially 
          ! ... filled by the topography
          flx = flx * ivf
!
          rhg(ijk) = - flx
!
! ... Same procedure carried out for solids
!
          DO is=1, nsolid
           IF (rlk(imjk,is) * inrl(is) <= 1.D-9) THEN
             esfx = b_e * esfe(ijk, is)
           ELSE
             esfx = b_e * esfe(ijk, is) - b_w * esfe(imjk, is)
           END IF

           IF (rlk(ijkm,is) * inrl(is) <= 1.D-9) THEN
             esfz = b_t * esft(ijk, is)
           ELSE
             esfz = b_t * esft(ijk, is) - b_b * esft(ijkm, is)
           END IF

           IF (job_type == JOB_TYPE_3D) THEN
             IF (rlk(ijmk,is) * inrl(is) <= 1.D-9) THEN
               esfy = b_n * esfn(ijk, is)
             ELSE
               esfy = b_n * esfn(ijk, is) - b_s * esfn(ijmk, is)
             END IF
           END IF
!
            flx = dt * indx(i) * esfx * inr(i) +  &
                  dt * indy(j) * esfy          +  &
                  dt * indz(k) * esfz
            !
            ! ... volume correction for cells partially 
            ! ... filled by the topography
            flx = flx * ivf
!
            rhs(ijk,is) = - flx + sink(ijk,is)*dt*ivf*sies(ijk,is)
            !
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
      IF (job_type == JOB_TYPE_3D) THEN
        DEALLOCATE(egfn)
        DEALLOCATE(esfn)
        DEALLOCATE(hgfn)
        DEALLOCATE(hsfn)
      END IF
!
      RETURN
      END SUBROUTINE htilde
!----------------------------------------------------------------------
      SUBROUTINE compute_all_fluxes_htilde

      USE atmospheric_conditions, ONLY: gravx, gravy, gravz
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE convective_fluxes_sc, ONLY: fsc, muscl_fsc, ctu1_fsc, ctu2_fsc, ctu3_fsc
      USE diffusive_fluxes, ONLY: hotc
      USE dimensions, ONLY: nsolid
      USE domain_mapping, ONLY: ncint, ncdom
      USE domain_mapping, ONLY: data_exchange
      USE flux_limiters, ONLY: muscl, ctu
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: sieg, sies, tg, ts
      USE gas_solid_temperature, ONLY: siegn
      USE gas_solid_viscosity, ONLY: kapg
      USE gas_solid_viscosity, ONLY: gas_viscosity, part_viscosity
      USE grid, ONLY: flag, fluid
      USE particles_constants, ONLY: inrl, kap
      USE pressure_epsilon, ONLY: ep
      USE set_indexes, ONLY: stencil, cte, subscr
      USE set_indexes, ONLY: ctu1_subscr, ctu2_subscr, ctu3_subscr
      USE set_indexes, ONLY: first_nb, first_rnb, third_nb, third_rnb
      USE set_indexes, ONLY: ctu1_nb, ctu1_rnb, ctu2_nb, ctu2_rnb, ctu3_nb, ctu3_rnb
      USE set_indexes, ONLY: imjk, ijmk, ijkm, OPERATOR( * )
      USE turbulence_model, ONLY: kapgt, iturb

      IMPLICIT NONE
!
      REAL*8 :: zero
      INTEGER :: is
      INTEGER :: ijk
      TYPE(stencil) :: u, v, w, dens, enth
      TYPE(stencil) :: eps, temp, kappa
      LOGICAL :: compute
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
        compute  = BTEST(flag(ijk),0)
        IF( compute ) THEN
          CALL subscr(ijk)
          IF (ctu > 0) CALL ctu1_subscr(ijk)
          IF (ctu > 1) CALL ctu2_subscr(ijk)
          IF (ctu > 2) CALL ctu3_subscr(ijk)
!
! ... Here compute convective and diffusive fluxes (gas)
!
          IF (job_type == JOB_TYPE_2D) THEN
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
            IF (muscl > 0 .AND. flag(ijk)==fluid) THEN
              CALL third_nb(enth,sieg,ijk)
              CALL third_nb(dens,rgp,ijk)
              CALL muscl_fsc(egfe(ijk), egft(ijk),    & 
                             enth, dens, u, w, ijk)
            END IF

            IF (ctu > 0) THEN
!
! ... First order Corner Transport Upwind correction (step 2)
!
              CALL ctu1_nb(enth,sieg,ijk)
              CALL ctu1_nb(dens,rgp,ijk)
              CALL ctu1_rnb(u,ug,ijk)
              CALL ctu1_rnb(w,wg,ijk)
              CALL ctu1_fsc(egfe(ijk), egft(ijk),    &
                      dens, enth, u, w, ijk)
!
! ... Second order Corner Transport Upwind correction (step 3)
!
              IF (ctu > 1 .AND. muscl == 0 .AND. flag(ijk)==fluid) THEN
                CALL ctu2_nb(enth,sieg,ijk)
                CALL ctu2_nb(dens,rgp,ijk)
                CALL ctu2_rnb(u,ug,ijk)
                CALL ctu2_rnb(w,wg,ijk)
                CALL ctu2_fsc(egfe(ijk), egft(ijk),    &
                        dens, enth, u, w, ijk)
!
! ... Third order Corner Transport Upwind correction (step 4)
!
                IF (ctu > 2) THEN
                  CALL ctu3_nb(enth,sieg,ijk)
                  CALL ctu3_nb(dens,rgp,ijk)
                  CALL ctu3_rnb(u,ug,ijk)
                  CALL ctu3_rnb(w,wg,ijk)
                  CALL ctu3_fsc(egfe(ijk), egft(ijk),    &
                          dens, enth, u, w, ijk)
                END IF
              END IF
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

          ELSE IF (job_type == JOB_TYPE_3D) THEN

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
            IF (muscl > 0 .AND. flag(ijk)==fluid) THEN
              CALL third_nb(enth,sieg,ijk)
              CALL third_nb(dens,rgp,ijk)
              CALL muscl_fsc(egfe(ijk), egfn(ijk), egft(ijk),    &
                             enth, dens, u, v, w, ijk)
            END IF
!
            IF (ctu > 0) THEN
!
! ... First order Corner Transport Upwind correction
!
              CALL ctu1_nb(enth,sieg,ijk)
              CALL ctu1_nb(dens,rgp,ijk)
              CALL ctu1_rnb(u,ug,ijk)
              CALL ctu1_rnb(v,ug,ijk)
              CALL ctu1_rnb(w,wg,ijk)
              CALL ctu1_fsc(egfe(ijk), egfn(ijk), egft(ijk),    &
                      enth, dens, u, v, w, ijk)
!
! ... Second order Corner Transport Upwind correction
!
              IF (ctu > 1 .AND. muscl == 0 .AND. flag(ijk) == fluid) THEN
                CALL ctu2_nb(enth,sieg,ijk)
                CALL ctu2_nb(dens,rgp,ijk)
                CALL ctu2_rnb(u,ug,ijk)
                CALL ctu2_rnb(v,ug,ijk)
                CALL ctu2_rnb(w,wg,ijk)
                CALL ctu2_fsc(egfe(ijk), egfn(ijk), egft(ijk),    &
                        enth, dens, u, v, w, ijk)
!
! ... Third order Corner Transport Upwind correction
!
                IF (ctu > 2) THEN
                  CALL ctu3_nb(enth,sieg,ijk)
                  CALL ctu3_nb(dens,rgp,ijk)
                  CALL ctu3_rnb(u,ug,ijk)
                  CALL ctu3_rnb(v,ug,ijk)
                  CALL ctu3_rnb(w,wg,ijk)
                  CALL ctu3_fsc(egfe(ijk), egfn(ijk), egft(ijk),    &
                          enth, dens, u, v, w, ijk)
                END IF
              END IF
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
          egfe(ijk) = egfe(ijk) - hgfe(ijk) + rgp(ijk)*gravx*ug(ijk)
          egft(ijk) = egft(ijk) - hgft(ijk) + rgp(ijk)*gravz*wg(ijk)
          IF ( .NOT.BTEST(flag(imjk),0) ) egfe(imjk) = egfe(imjk) - hgfe(imjk)
          IF ( .NOT.BTEST(flag(ijkm),0) ) egft(ijkm) = egft(ijkm) - hgft(ijkm)
!
          IF (job_type == JOB_TYPE_3D) THEN
            egfn(ijk) = egfn(ijk) - hgfn(ijk) + rgp(ijk)*gravy*vg(ijk)
            IF ( .NOT.BTEST(flag(ijmk),0) ) egfn(ijmk) = egfn(ijmk) - hgfn(ijmk)
          END IF
!
! ... Here compute convective and diffusive fluxes (particles)
!
          DO is=1, nsolid
!
            IF (job_type == JOB_TYPE_2D) THEN

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
              IF (muscl > 0 .AND. flag(ijk)==fluid) THEN
                CALL third_nb(enth,sies(:,is),ijk)
                CALL third_nb(dens,rlk(:,is),ijk)
                CALL muscl_fsc(esfe(ijk, is), esft(ijk, is),  &
                               enth, dens, u, w, ijk)
              END IF

              IF (ctu > 0) THEN
!
! ... First order Corner Transport Upwind correction (step 2)
!
                CALL ctu1_nb(enth,sies(:,is),ijk)
                CALL ctu1_nb(dens,rlk(:,is),ijk)
                CALL ctu1_rnb(u, us(:,is),ijk)
                CALL ctu1_rnb(w, ws(:,is),ijk)
                CALL ctu1_fsc(esfe(ijk, is), esft(ijk, is),  &
                        dens, enth, u, w, ijk)
!
! ... Second order Corner Transport Upwind correction (step 3)
!
                IF (ctu > 1 .AND. muscl == 0 .AND. flag(ijk) ==fluid) THEN
                  CALL ctu2_nb(enth,sies(:,is),ijk)
                  CALL ctu2_nb(dens,rlk(:,is),ijk)
                  CALL ctu2_rnb(u, us(:,is),ijk)
                  CALL ctu2_rnb(w, ws(:,is),ijk)
                  CALL ctu2_fsc(esfe(ijk, is), esft(ijk, is),  &
                          dens, enth, u, w, ijk)
!
! ... Third order Corner Transport Upwind (step 4)
!
                  IF (ctu > 2) THEN
                    CALL ctu3_nb(enth,sies(:,is),ijk)
                    CALL ctu3_nb(dens,rlk(:,is),ijk)
                    CALL ctu3_rnb(u, us(:,is),ijk)
                    CALL ctu3_rnb(w, ws(:,is),ijk)
                    CALL ctu3_fsc(esfe(ijk, is), esft(ijk, is),  &
                            dens, enth, u, w, ijk)
                  END IF
                END IF
              END IF
!
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

            ELSE IF (job_type == JOB_TYPE_3D) THEN

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
              IF (muscl > 0 .AND. flag(ijk)==fluid) THEN
                CALL third_nb(enth,sies(:,is),ijk)
                CALL third_nb(dens,rlk(:,is),ijk)
                CALL muscl_fsc(esfe(ijk, is), esfn(ijk, is), esft(ijk, is),  &
                               enth, dens, u, v, w, ijk)
              END IF
!
              IF (ctu > 0) THEN
!
! ... First order Corner Transport Upwind
!
                CALL ctu1_nb(enth,sies(:,is),ijk)
                CALL ctu1_nb(dens,rlk(:,is),ijk)
                CALL ctu1_rnb(u, us(:,is),ijk)
                CALL ctu1_rnb(v, vs(:,is),ijk)
                CALL ctu1_rnb(w, ws(:,is),ijk)
                CALL ctu1_fsc(esfe(ijk, is), esfn(ijk, is), esft(ijk, is),  &
                        enth, dens, u, v, w, ijk)
!
! ... Second order Corner Transport Upwind
!
                IF (ctu > 1) THEN
                  CALL ctu2_nb(enth,sies(:,is),ijk)
                  CALL ctu2_nb(dens,rlk(:,is),ijk)
                  CALL ctu2_rnb(u, us(:,is),ijk)
                  CALL ctu2_rnb(v, vs(:,is),ijk)
                  CALL ctu2_rnb(w, ws(:,is),ijk)
                  CALL ctu2_fsc(esfe(ijk, is), esfn(ijk, is), esft(ijk, is),  &
                          enth, dens, u, v, w, ijk)
!
! ... Third order Corner Transport Upwind
!
                  IF (ctu > 2) THEN
                    CALL ctu3_nb(enth,sies(:,is),ijk)
                    CALL ctu3_nb(dens,rlk(:,is),ijk)
                    CALL ctu3_rnb(u, us(:,is),ijk)
                    CALL ctu3_rnb(v, vs(:,is),ijk)
                    CALL ctu3_rnb(w, ws(:,is),ijk)
                    CALL ctu3_fsc(esfe(ijk, is), esfn(ijk, is), esft(ijk, is),  &
                            enth, dens, u, v, w, ijk)
                  END IF
                END IF
              END IF
!
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
            esfe(ijk, is) = esfe(ijk, is) - hsfe(ijk,is) + rlk(ijk,is)*gravx*us(ijk,is)
            esft(ijk, is) = esft(ijk, is) - hsft(ijk, is) + rlk(ijk,is)*gravz*ws(ijk,is)
            IF ( .NOT.BTEST(flag(imjk),0) ) esfe(imjk,is) = esfe(imjk,is) - hsfe(imjk,is)
            IF ( .NOT.BTEST(flag(ijkm),0) ) esft(ijkm,is) = esft(ijkm,is) - hsft(ijkm,is)
!
            IF (job_type == JOB_TYPE_3D) THEN
              esfn(ijk, is) = esfn(ijk, is) - hsfn(ijk, is) + rlk(ijk,is)*gravy*vs(ijk,is)
              IF ( .NOT.BTEST(flag(ijmk),0) ) esfn(ijmk,is)=esfn(ijmk,is)-hsfn(ijmk,is)
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

      IF (job_type == JOB_TYPE_3D) THEN
        CALL data_exchange(egfn)
        CALL data_exchange(esfn)
      END IF

      RETURN
      END SUBROUTINE compute_all_fluxes_htilde
!----------------------------------------------------------------------
      SUBROUTINE test_fluxes_htilde
      USE domain_mapping, ONLY: ncint, meshinds
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm
      USE gas_solid_velocity, ONLY: ug
      USE io_files, ONLY: tempunit

      IMPLICIT NONE
      INTEGER :: ijk, imesh, i, j, k

      OPEN(UNIT=tempunit,FILE='pdac.fl',STATUS='UNKNOWN')
      WRITE(tempunit,*) 'Test enthalpy fluxes ...'
      WRITE(tempunit,*)
      DO ijk = 1, ncint
        CALL subscr(ijk)
        CALL meshinds(ijk,imesh,i,j,k)

        !WRITE(tempunit,100) ijk, i, j, k, egfe(ijk), egfn(ijk), egft(ijk)
        WRITE(tempunit,101) egfe(ijk), egfn(ijk), egft(ijk)

      END DO

 100  FORMAT(4(I5),3(G20.10E2))
 101  FORMAT(3(G20.10E2))
      CLOSE(tempunit)
      
      RETURN
      END SUBROUTINE test_fluxes_htilde
!----------------------------------------------------------------------
      END MODULE tilde_energy
!-----------------------------------------------------------
