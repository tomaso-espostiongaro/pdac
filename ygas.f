!----------------------------------------------------------------------
      MODULE gas_components
!----------------------------------------------------------------------
      IMPLICIT NONE
!
! ... convective concentration fluxes 
!
      REAL*8, ALLOCATABLE :: yfe(:,:), yfn(:,:), yft(:,:)
      
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE ygas
!
      USE control_flags, ONLY: job_type, lpr
      USE dimensions
      USE domain_decomposition, ONLY: ncint, ncdom, meshinds
      USE eos_gas, ONLY: rgpgc, rgpgcn, ygc, xgc, mole
      USE gas_constants, ONLY: default_gas, gas_type
      USE grid, ONLY: dx, dy, dz, indx, indy, indz, inx
      USE grid, ONLY: fl_l
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm
      USE time_parameters, ONLY: dt,time
!
      IMPLICIT NONE
!
      REAL*8 :: yfw, yfs, yfb, yfx, yfy, yfz
      REAL*8 :: ygcdfg
      REAL*8 :: rgp, rgpinv, rgpgc_tmp
      INTEGER :: i, j, k, imesh
      INTEGER :: ijk
      INTEGER :: ig, dfg
!
      ALLOCATE(yfe(ncdom,ngas), yft(ncdom,ngas))
      yfe = 0.D0; yft = 0.D0

      IF (job_type == '3D') THEN
        ALLOCATE(yfn(ncdom,ngas))
        yfn = 0.D0
      END IF
! 
! ... Convective mass fluxes
!
      CALL compute_all_fluxes

      DO ijk = 1, ncint
       IF( fl_l(ijk) == 1 ) THEN
         CALL meshinds(ijk,imesh,i,j,k)
         CALL subscr(ijk)
	   
         yfx = 0.D0
         yfy = 0.D0
         yfz = 0.D0
         rgp = 0.D0

         DO ig=1,ngas
	   
	   yfw = yfe(imjk,ig)
	   yfb = yft(ijkm,ig)

           yfx = yfe(ijk,ig) - yfw
	   yfz = yft(ijk,ig) - yfb

           IF (job_type == '3D') THEN
  	     yfs = yfn(ijmk,ig)
  	     yfy = yfn(ijk,ig) - yfs
           END IF

           rgpgc_tmp = rgpgcn(ijk,ig)
	   rgpgc_tmp = rgpgc_tmp - dt * indx(i) * yfx * inx(i)     
           rgpgc_tmp = rgpgc_tmp - dt * indy(j) * yfy              
	   rgpgc_tmp = rgpgc_tmp - dt * indz(k) * yfz
 
           IF (lpr > 1) THEN
             IF(rgpgc(ijk,ig) < 0.D0) THEN
               WRITE(8,*) 'Warning!: gas mass is not conserved'
               WRITE(8,128) time, ijk, ig 
               WRITE(8,*) rgpgcn(ijk,ig),  rgpgc(ijk,ig)
 128           FORMAT('Time= ',F8.3,' Cell= ',I6, ' Specie= ',I2)
               rgpgc(ijk,ig) = 0.D0
             END IF
           END IF

           rgpgc(ijk,ig) = MAX(rgpgc_tmp,0.D0)
           
           rgp = rgp + rgpgc(ijk,ig)

         END DO

         rgpinv = 1.0D0/rgp 
!
! ... Update the mass fractions (with the closure relation constraint)
!
         ygcdfg = 1.D0
         DO ig=1,ngas
           IF (gas_type(ig) /= default_gas) THEN
             ygc(ig,ijk) = rgpgc(ijk,ig) * rgpinv
             ygcdfg = ygcdfg - ygc(ig,ijk)
           ELSE 
             dfg = ig
           END IF
         END DO
         ygc(dfg,ijk) = ygcdfg
!
! ... Update molar fractions
!
         CALL mole( xgc(:,ijk), ygc(:,ijk) )
!
       END IF
      END DO
!
      DEALLOCATE(yfe, yft)
      IF (job_type == '3D') DEALLOCATE(yfn)
!
      RETURN
      END SUBROUTINE ygas
!----------------------------------------------------------------------
      SUBROUTINE compute_all_fluxes
!
      USE control_flags, ONLY: job_type
      USE convective_fluxes_sc, ONLY: fsc, muscl_fsc
      USE dimensions, ONLY: ngas
      USE domain_decomposition, ONLY: ncint, data_exchange
      USE eos_gas, ONLY: rgpgc
      USE flux_limiters, ONLY: muscl
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE grid, ONLY: fl_l
      USE set_indexes, ONLY: stencil, cte
      USE set_indexes, ONLY: first_nb, first_rnb, third_nb, third_rnb
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm

      IMPLICIT NONE

      INTEGER :: ijk, ig
      TYPE(stencil) :: one, field, u, v, w
!
      CALL data_exchange(rgpgc)

      one  = cte(1.D0)
!
      DO ijk = 1, ncint
       IF( fl_l(ijk) == 1 ) THEN
         CALL subscr(ijk)

         DO ig=1,ngas
           
           IF (job_type == '2D') THEN
!
! ... Compute fluxes by using First Order Upwinds
!
             ! ... Assemble the first order computational stencils
	     CALL first_nb(field,rgpgc(:,ig),ijk)
             CALL first_rnb(u,ug,ijk)
	     CALL first_rnb(w,wg,ijk)

  	     CALL fsc(yfe(ijk,ig), yft(ijk,ig),      &
                      yfe(imjk,ig), yft(ijkm,ig),   &
                      one, field, u, w, ijk)
!
! ... Second order MUSCL correction
!
             IF (muscl > 0) THEN

               CALL third_nb(field,rgpgc(:,ig),ijk)
  	       CALL muscl_fsc(yfe(ijk,ig), yft(ijk,ig),      &
                            one, field, u, w, ijk)

             END IF

           ELSE IF (job_type == '3D') THEN
!
! ... Compute fluxes by using First Order Upwinds
!
             ! ... Assemble the first order computational stencils
	     CALL first_nb(field,rgpgc(:,ig),ijk)
             CALL first_rnb(u,ug,ijk)
	     CALL first_rnb(v,vg,ijk)
	     CALL first_rnb(w,wg,ijk)

  	     CALL fsc(yfe(ijk,ig), yfn(ijk,ig), yft(ijk,ig),      &
                      yfe(imjk,ig), yfn(ijmk,ig), yft(ijkm,ig),   &
                      one, field, u, v, w, ijk)
!
! ... Second order MUSCL correction
!
             IF (muscl > 0) THEN

               CALL third_nb(field,rgpgc(:,ig),ijk)
  	       CALL muscl_fsc(yfe(ijk,ig), yfn(ijk,ig), yft(ijk,ig),      &
                              one, field, u, v, w, ijk)

             END IF

           END IF
         END DO
         
       END IF
      END DO
      
      CALL data_exchange(yfe)
      CALL data_exchange(yft)
      IF (job_type == '3D') CALL data_exchange(yfn)

      END SUBROUTINE compute_all_fluxes
!----------------------------------------------------------------------
      END MODULE gas_components
!----------------------------------------------------------------------
