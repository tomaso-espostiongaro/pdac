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
      USE domain_decomposition, ONLY: myijk, data_exchange,myinds
      USE domain_decomposition, ONLY: ncint, ncdom, meshinds
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: sieg, sies, tg, ts
      USE gas_solid_viscosity, ONLY: kapg
      USE grid, ONLY: dx, dy, dz, indx, indy, indz, inx
      USE grid, ONLY: fl_l
      USE indijk_module
      USE particles_constants, ONLY: inrl
      USE pressure_epsilon, ONLY: ep
      USE set_indexes, ONLY: subscr, subscr_red, imjk, ijmk, ijkm
      USE set_indexes, ONLY: ijke, ijkw, ijkn, ijks, ijkt, ijkb
      USE time_parameters, ONLY: dt
      USE turbulence_model, ONLY: kapgt, iturb
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
! ... Compute East, North, and Top fluxes in every cell
! ... within the computational domain
! 
      IF ( job_type == '2D' ) THEN
        CALL compute_all_fluxes
      ELSE IF ( job_type == '3D' ) THEN
        CALL compute_all_fluxes_3d
      END IF
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
          !CALL subscr(ijk)
          !CALL subscr_red(ijk)
          IF( job_type == '2D' ) THEN

                 ijkm  = myijk( ip0_jp0_km1_, ijk )
                 imjk  = myijk( im1_jp0_kp0_, ijk )

                 ijke  = myinds(ip1_jp0_kp0_, ijk )
                 ijkt  = myinds(ip0_jp0_kp1_, ijk )
                 ijkw  = myinds(im1_jp0_kp0_, ijk )
                 ijkb  = myinds(ip0_jp0_km1_, ijk )

          ELSE IF( job_type == '3D' ) THEN

                 imjk   = myijk( im1_jp0_kp0_ , ijk )
                 ijmk   = myijk( ip0_jm1_kp0_ , ijk )
                 ijkm   = myijk( ip0_jp0_km1_ , ijk )

                 ijke = myinds( ip1_jp0_kp0_ , ijk )
                 ijkw = myinds( im1_jp0_kp0_ , ijk )
                 ijkn = myinds( ip0_jp1_kp0_ , ijk )
                 ijks = myinds( ip0_jm1_kp0_ , ijk )
                 ijkt = myinds( ip0_jp0_kp1_ , ijk )
                 ijkb = myinds( ip0_jp0_km1_ , ijk )

          END IF

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
      SUBROUTINE compute_all_fluxes

      USE control_flags, ONLY: job_type
      USE convective_fluxes_sc, ONLY: fsc
      USE diffusive_fluxes, ONLY: hotc
      USE dimensions, ONLY: nsolid
      USE domain_decomposition, ONLY: ncint, ncdom
      USE domain_decomposition, ONLY: data_exchange
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: sieg, sies, tg, ts
      USE gas_solid_viscosity, ONLY: kapg
      USE gas_solid_viscosity, ONLY: gas_viscosity, part_viscosity
      USE grid, ONLY: fl_l
      USE particles_constants, ONLY: inrl, kap
      USE pressure_epsilon, ONLY: ep
      USE set_indexes, ONLY: nb, rnb, stencil, cte
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes
      USE turbulence_model, ONLY: kapgt

      IMPLICIT NONE
!
      REAL*8 :: zero
      INTEGER :: is
      INTEGER :: ijk
      TYPE(stencil) :: u, v, w, dens, enth
      TYPE(stencil) :: eps, temp, kappa

      DO ijk = 1, ncint
        IF(fl_l(ijk) == 1) THEN
          !CALL subscr(ijk)
          CALL subscr_red(ijk)
!
! ... compute convective and diffusive fluxes (gas)
!
! ... assemble the computational stencils
          CALL nb_13(enth,sieg,ijk)
          CALL nb_13(dens,rgp,ijk)
 
          !CALL nb(kappa,kapgt,ijk)
          kappa%c = kapgt ( ijk )
          kappa%w = kapgt ( ijkw )
          kappa%s = kapgt ( ijks )
          kappa%b = kapgt ( ijkb )
          kappa%e = kapgt ( ijke )
          kappa%n = kapgt ( ijkn )
          kappa%t = kapgt ( ijkt )
          !CALL nb(eps,ep,ijk)
          eps%c = ep ( ijk )
          eps%w = ep ( ijkw )
          eps%s = ep ( ijks )
          eps%b = ep ( ijkb )
          eps%e = ep ( ijke )
          eps%n = ep ( ijkn )
          eps%t = ep ( ijkt )
          !CALL nb(temp,tg,ijk)
          temp%c = tg ( ijk )
          temp%w = tg ( ijkw )
          temp%s = tg ( ijks )
          temp%b = tg ( ijkb )
          temp%e = tg ( ijke )
          temp%n = tg ( ijkn )
          temp%t = tg ( ijkt )


          !CALL rnb(u,ug,ijk)
          !CALL rnb(w,wg,ijk)
          u%w = ug ( imjk )
          u%c = ug ( ijk )
          w%b = wg ( ijkm )
          w%c = wg ( ijk )  
      
          IF (job_type == '2D') THEN

            CALL fsc(egfe(ijk), egft(ijk),    &
                     egfe(imjk), egft(ijkm),    &
                     enth, dens, u, w, ijk)
!
            IF (gas_viscosity) THEN
              CALL hotc(hgfe(ijk), zero, hgft(ijk),        &
                        hgfe(imjk), zero, hgft(ijkm),      &
                        eps, temp, kappa, ijk)
            END IF

          ELSE IF (job_type == '3D') THEN

            !CALL rnb(v,vg,ijk)
            v%s = vg ( ijmk )
            v%c = vg ( ijk )

            CALL fsc(egfe(ijk), egfn(ijk), egft(ijk),    &
                     egfe(imjk), egfn(ijmk), egft(ijkm),    &
                     enth, dens, u, v, w, ijk)
!
            IF (gas_viscosity) THEN
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
          DO is=1, nsolid
!
! ... convective and diffusive fluxes (particles)
!
! ... assemble computational stencils
            CALL nb(enth,sies(:,is),ijk)
            CALL nb(dens,rlk(:,is),ijk)
            CALL nb(temp,ts(:,is),ijk)
            !CALL rnb(u, us(:,is),ijk)
            !CALL rnb(w, ws(:,is),ijk)
            u%w = us ( imjk , is)
            u%c = us ( ijk , is)
            w%b = ws ( ijkm , is)
            w%c = ws ( ijk , is ) 
            eps   = inrl(is) * dens
            kappa = cte(kap(is))
            
            IF (job_type == '2D') THEN

              CALL fsc(esfe(ijk, is), esft(ijk, is),  &
                       esfe(imjk, is), esft(ijkm, is),  &
                       enth, dens, u, w, ijk)
!
              IF (part_viscosity) THEN
                CALL hotc(hsfe(ijk,is), zero, hsft(ijk, is),    &
                          hsfe(imjk,is), zero, hsft(ijkm, is),    &
                          eps, temp, kappa, ijk)
              END IF

            ELSE IF (job_type == '3D') THEN

              !CALL rnb(v, vs(:,is),ijk)
              v%s = vs ( ijmk , is )
              v%c = vs ( ijk , is )
              CALL fsc(esfe(ijk, is), esfn(ijk, is), esft(ijk, is),  &
                       esfe(imjk, is), esfn(ijmk, is), esft(ijkm, is),  &
                       enth, dens, u, v, w, ijk)
!
              IF (part_viscosity) THEN
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
      SUBROUTINE compute_all_fluxes_3d

      USE convective_fluxes_sc, ONLY: fsc, fsc_1st
      USE diffusive_fluxes, ONLY: hotc
      USE dimensions, ONLY: nsolid
      USE domain_decomposition, ONLY: ncint, ncdom
      USE domain_decomposition, ONLY: data_exchange, myinds
      USE flux_limiters, ONLY: muscl
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: sieg, sies, tg, ts
      USE gas_solid_viscosity, ONLY: kapg
      USE gas_solid_viscosity, ONLY: gas_viscosity, part_viscosity
      USE grid, ONLY: fl_l
      USE indijk_module
      USE particles_constants, ONLY: inrl, kap
      USE pressure_epsilon, ONLY: ep
      USE set_indexes, ONLY: nb, rnb, stencil, cte
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes
      USE turbulence_model, ONLY: kapgt

      IMPLICIT NONE
!
      REAL*8 :: zero
      INTEGER :: is
      INTEGER :: ijk
      TYPE(stencil) :: u, v, w, dens, enth
      TYPE(stencil) :: eps, temp, kappa

      DO ijk = 1, ncint
        IF(fl_l(ijk) == 1) THEN
          !CALL subscr(ijk)
          !CALL subscr_red(ijk)

             
                 imjk   = myijk( im1_jp0_kp0_ , ijk )
                 ijmk   = myijk( ip0_jm1_kp0_ , ijk )
                 ijkm   = myijk( ip0_jp0_km1_ , ijk )

                 ijke = myinds( ip1_jp0_kp0_ , ijk )
                 ijkw = myinds( im1_jp0_kp0_ , ijk )
                 ijkn = myinds( ip0_jp1_kp0_ , ijk )
                 ijks = myinds( ip0_jm1_kp0_ , ijk )
                 ijkt = myinds( ip0_jp0_kp1_ , ijk )
                 ijkb = myinds( ip0_jp0_km1_ , ijk )

             

!
! ... compute convective and diffusive fluxes (gas)
!
! ... assemble the computational stencils
          CALL nb_13(enth,sieg,ijk)
          CALL nb_13(dens,rgp,ijk)
          !CALL nb(kappa,kapgt,ijk)
          kappa%c = kapgt ( ijk )
          kappa%w = kapgt ( ijkw )
          kappa%s = kapgt ( ijks )
          kappa%b = kapgt ( ijkb )
          kappa%e = kapgt ( ijke )
          kappa%n = kapgt ( ijkn )
          kappa%t = kapgt ( ijkt )
          !CALL nb(eps,ep,ijk)
          eps%c = ep ( ijk )
          eps%w = ep ( ijkw )
          eps%s = ep ( ijks )
          eps%b = ep ( ijkb )
          eps%e = ep ( ijke )
          eps%n = ep ( ijkn )
          eps%t = ep ( ijkt )
          !CALL nb(temp,tg,ijk)
          temp%c = tg ( ijk )
          temp%w = tg ( ijkw )
          temp%s = tg ( ijks )
          temp%b = tg ( ijkb )
          temp%e = tg ( ijke )
          temp%n = tg ( ijkn )
          temp%t = tg ( ijkt )

          !CALL rnb(u,ug,ijk)
          !CALL rnb(w,wg,ijk)
          u%w = ug ( imjk )
          u%c = ug ( ijk )
          w%b = wg ( ijkm )
          w%c = wg ( ijk )  
          !CALL rnb(v,vg,ijk)
          v%s = vg ( ijmk )
          v%c = vg ( ijk )

            IF ( muscl == 0 ) THEN
              CALL fsc_1st(egfe(ijk), egfn(ijk), egft(ijk),    &
                     egfe(imjk), egfn(ijmk), egft(ijkm),    &
                     enth, dens, u, v, w)
            ELSE 
              CALL fsc(egfe(ijk), egfn(ijk), egft(ijk),    &
                     egfe(imjk), egfn(ijmk), egft(ijkm),    &
                     enth, dens, u, v, w, ijk)
            END IF
!
            IF (gas_viscosity) THEN
              CALL hotc(hgfe(ijk), hgfn(ijk), hgft(ijk),         &
                        hgfe(imjk), hgfn(ijmk), hgft(ijkm),      &
                        eps, temp, kappa, ijk)
            END IF

!
            IF (gas_viscosity) THEN
              egfe(ijk) = egfe(ijk) - hgfe(ijk)
              egft(ijk) = egft(ijk) - hgft(ijk)
!
              IF (fl_l(imjk) /= 1) egfe(imjk) = egfe(imjk) - hgfe(imjk)
              IF (fl_l(ijkm) /= 1) egft(ijkm) = egft(ijkm) - hgft(ijkm)         
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
            !CALL rnb(u, us(:,is),ijk)
            !CALL rnb(v, vs(:,is),ijk)
            !CALL rnb(w, ws(:,is),ijk)           
            u%w = us ( imjk , is)
            u%c = us ( ijk , is)
            v%s = vs ( ijmk , is )
            v%c = vs ( ijk , is )
            w%b = ws ( ijkm , is)
            w%c = ws ( ijk , is ) 
            eps   = inrl(is) * dens
            kappa = cte(kap(is))
            
              IF ( muscl == 0 ) THEN
                CALL fsc_1st(esfe(ijk, is), esfn(ijk, is), esft(ijk, is),  &
                       esfe(imjk, is), esfn(ijmk, is), esft(ijkm, is),  &
                       enth, dens, u, v, w)
              ELSE
                CALL fsc(esfe(ijk, is), esfn(ijk, is), esft(ijk, is),  &
                       esfe(imjk, is), esfn(ijmk, is), esft(ijkm, is),  &
                       enth, dens, u, v, w, ijk)
              END IF
!
              IF (part_viscosity) THEN
                CALL hotc(hsfe(ijk,is), hsfn(ijk, is), hsft(ijk, is),    &
                          hsfe(imjk,is), hsfn(ijmk, is), hsft(ijkm, is),    &
                          eps, temp, kappa, ijk)
              END IF
!
              IF (part_viscosity) THEN
                esfe(ijk, is) = esfe(ijk, is) - hsfe(ijk,is)
                esft(ijk, is) = esft(ijk, is) - hsft(ijk, is)
                IF (fl_l(imjk) /= 1) esfe(imjk,is) = esfe(imjk,is) - hsfe(imjk,is)
                IF (fl_l(ijkm) /= 1) esft(ijkm,is) = esft(ijkm,is) - hsft(ijkm,is)             
                esfn(ijk, is) = esfn(ijk, is) - hsfn(ijk, is)
                IF (fl_l(ijmk)/=1) esfn(ijmk,is)=esfn(ijmk,is)-hsfn(ijmk,is)
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
      CALL data_exchange(egfn)
      CALL data_exchange(esfn)

      RETURN
      END SUBROUTINE compute_all_fluxes_3d
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
