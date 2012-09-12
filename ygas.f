!----------------------------------------------------------------------
      MODULE gas_components
!----------------------------------------------------------------------
      USE io_files, ONLY: testunit
      IMPLICIT NONE
!
! ... convective concentration fluxes 
!
      REAL*8, ALLOCATABLE :: yfe(:,:), yfn(:,:), yft(:,:)
      REAL*8, ALLOCATABLE :: rgpgcn(:,:)
      
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_species
      USE dimensions
      USE domain_mapping, ONLY: ncint
      IMPLICIT NONE
      
      ALLOCATE(rgpgcn(ncint,ngas))
      rgpgcn = 0.D0

      RETURN
      END SUBROUTINE allocate_species
!----------------------------------------------------------------------
      SUBROUTINE ygas
!
      USE control_flags, ONLY: job_type, lpr
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions
      USE domain_mapping, ONLY: ncint, ncdom, meshinds
      USE eos_gas, ONLY: ygc, xgc, mole
      USE gas_constants, ONLY: gas_type
      USE gas_solid_density, ONLY: rgp, rgpn
      USE grid, ONLY: dx, dy, dz, indx, indy, indz, inr
      USE grid, ONLY: flag
      USE immersed_boundaries, ONLY: faces, immb
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm
      USE time_parameters, ONLY: dt,time, sweep
!
      IMPLICIT NONE
!
      REAL*8 :: yfw, yfs, yfb, yfx, yfy, yfz
      REAL*8 :: ygcdfg
      REAL*8 :: rgpt, rgpinv
      INTEGER :: i, j, k, imesh
      INTEGER :: ijk
      INTEGER :: ig, dfg
      INTEGER :: b_e, b_w, b_t, b_b, b_n, b_s
      REAL*8 :: ivf

      REAL*8, ALLOCATABLE :: rgpgc(:)
      LOGICAL :: compute
!
      IF (ngas == 1) THEN
        ygc(:,1) = 1.D0
        xgc(:,1) = 1.D0
        RETURN
      END IF

      ALLOCATE(rgpgc(ngas)); rgpgc = 0.D0
!
      ALLOCATE(yfe(ncdom,ngas), yft(ncdom,ngas))
      yfe = 0.D0; yft = 0.D0

      IF (job_type == JOB_TYPE_3D) THEN
        ALLOCATE(yfn(ncdom,ngas))
        yfn = 0.D0
      END IF
! 
! ... Compute the convective mass fluxes
!
      CALL compute_all_fluxes
!
! ... Mass conservation of each gas component is solved
! ... explicitly in each cell
!
      DO ijk = 1, ncint
       compute = BTEST(flag(ijk),0)
       IF( compute ) THEN
         !
         ! ... Compute the volumes partially filled by the
         ! ... topography
         CALL faces(ijk, b_e, b_w, b_t, b_b, b_n, b_s, ivf)

         CALL meshinds(ijk,imesh,i,j,k)
         CALL subscr(ijk)

         yfx = 0.D0
         yfy = 0.D0
         yfz = 0.D0
         rgpt = 0.D0

         DO ig=1,ngas
	   
	   yfw = yfe(imjk,ig)
	   yfb = yft(ijkm,ig)

           yfx = b_e * yfe(ijk,ig) - b_w * yfw
	   yfz = b_t * yft(ijk,ig) - b_b * yfb

           IF (job_type == JOB_TYPE_3D) THEN
  	     yfs = yfn(ijmk,ig)
  	     yfy = b_n * yfn(ijk,ig) - b_s * yfs
           END IF

           rgpgc(ig) = rgpgcn(ijk,ig)
	   rgpgc(ig) = rgpgc(ig) - ivf * dt * indx(i) * yfx * inr(i)     
           rgpgc(ig) = rgpgc(ig) - ivf * dt * indy(j) * yfy              
	   rgpgc(ig) = rgpgc(ig) - ivf * dt * indz(k) * yfz

           IF(rgpgc(ig) < 0.D0) THEN
             IF (lpr > 2) THEN
               WRITE(testunit,128) ijk, ig 
               WRITE(testunit,*) rgpgc(ig)
 128           FORMAT(' Cell= ',I6, ' Specie= ',I2)
             END IF

             rgpgc(ig) = 0.D0

           END IF

           rgpt = rgpt + rgpgc(ig)

         END DO

         rgpinv = 1.0D0/rgpt 
!
! ... Update the mass fractions (with the closure relation constraint)
!
        ygcdfg = 1.D0
        DO ig = 1, ngas - 1
          ygc(ijk,ig) = rgpgc(ig) * rgpinv
          ygcdfg = ygcdfg - ygc(ijk,ig)
        END DO
        ygc(ijk,ngas) = ygcdfg
!
! ... Update molar fractions
!
        CALL mole( xgc(ijk,:), ygc(ijk,:) )
!
       END IF
      END DO
!
      DEALLOCATE(rgpgc)
      DEALLOCATE(yfe, yft)
      IF (job_type == JOB_TYPE_3D) DEALLOCATE(yfn)
!
      RETURN
      END SUBROUTINE ygas
!----------------------------------------------------------------------
      SUBROUTINE compute_all_fluxes
!
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE convective_fluxes_sc, ONLY: fsc, muscl_fsc
      USE dimensions, ONLY: ngas
      USE domain_mapping, ONLY: ncint, data_exchange
      USE eos_gas, ONLY: ygc
      USE flux_limiters, ONLY: muscl
      USE gas_solid_density, ONLY: rgp
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE grid, ONLY: flag, fluid
      USE set_indexes, ONLY: stencil
      USE set_indexes, ONLY: first_nb, first_rnb, third_nb, third_rnb
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm

      IMPLICIT NONE

      INTEGER :: ijk, ig
      TYPE(stencil) :: dens, conc, u, v, w
      LOGICAL :: compute
!
      CALL data_exchange(ygc)
!
      DO ijk = 1, ncint
        compute  = BTEST(flag(ijk),0)
       IF( compute ) THEN
         CALL subscr(ijk)

         DO ig=1,ngas
           
           IF (job_type == JOB_TYPE_2D ) THEN
!
! ... Compute fluxes by using First Order Upwinds
!
             ! ... Assemble the first order computational stencils
	     CALL first_nb(dens,rgp(:),ijk)
	     CALL first_nb(conc,ygc(:,ig),ijk)
             CALL first_rnb(u,ug,ijk)
	     CALL first_rnb(w,wg,ijk)

  	     CALL fsc(yfe(ijk,ig), yft(ijk,ig),      &
                      yfe(imjk,ig), yft(ijkm,ig),   &
                      dens, conc, u, w, ijk)
!
! ... Second order MUSCL correction
!
             IF (muscl > 0 .AND. flag(ijk)==fluid) THEN

               CALL third_nb(dens,rgp(:),ijk)
               CALL third_nb(conc,ygc(:,ig),ijk)
  	       CALL muscl_fsc(yfe(ijk,ig), yft(ijk,ig),      &
                            dens, conc, u, w, ijk)

             END IF

           ELSE IF (job_type == JOB_TYPE_3D) THEN
!
! ... Compute fluxes by using First Order Upwinds
!
             ! ... Assemble the first order computational stencils
	     CALL first_nb(dens,rgp(:),ijk)
	     CALL first_nb(conc,ygc(:,ig),ijk)
             CALL first_rnb(u,ug,ijk)
	     CALL first_rnb(v,vg,ijk)
	     CALL first_rnb(w,wg,ijk)

  	     CALL fsc(yfe(ijk,ig), yfn(ijk,ig), yft(ijk,ig),      &
                      yfe(imjk,ig), yfn(ijmk,ig), yft(ijkm,ig),   &
                      dens, conc, u, v, w, ijk)
!
! ... Second order MUSCL correction
!
             IF (muscl > 0 .AND. flag(ijk)==fluid) THEN

               CALL third_nb(dens,rgp(:),ijk)
               CALL third_nb(conc,ygc(:,ig),ijk)
  	       CALL muscl_fsc(yfe(ijk,ig), yfn(ijk,ig), yft(ijk,ig),      &
                              dens, conc, u, v, w, ijk)

             END IF

           END IF
         END DO
         
       END IF
      END DO
      
      CALL data_exchange(yfe)
      CALL data_exchange(yft)
      IF (job_type == JOB_TYPE_3D) CALL data_exchange(yfn)

      END SUBROUTINE compute_all_fluxes
!----------------------------------------------------------------------
      END MODULE gas_components
!----------------------------------------------------------------------
