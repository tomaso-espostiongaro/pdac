!----------------------------------------------------------------------
      MODULE gas_components
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
      REAL*8, ALLOCATABLE :: yfe(:,:), yfn(:,:), yft(:,:)
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE ygas
!
      USE convective_fluxes_sc, ONLY: fsc
      USE control_flags, ONLY: job_type
      USE dimensions
      USE eos_gas, ONLY: rgpgc, rgpgcn, ygc
      USE gas_constants, ONLY: default_gas
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE grid, ONLY: dx, dy, dz, indx, indy, indz
      USE grid, ONLY: ncint, ncdom, myijk, data_exchange
      USE grid, ONLY: fl_l, meshinds
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil, nb, rnb, cte
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm
      USE time_parameters, ONLY: dt,time
!
      IMPLICIT NONE
!
      TYPE(stencil) :: one, field, u, v, w
      REAL*8 :: yfw, yfs, yfb, yfx, yfy, yfz
      REAL*8 :: rgp
      INTEGER :: i, j, k, imesh
      INTEGER :: ijk
      INTEGER :: ig
!
      ALLOCATE(yfe(ncdom,ngas), yft(ncdom,ngas))
      yfe = 0.D0; yft = 0.D0

      IF (job_type == '3D') THEN
        ALLOCATE(yfn(ncdom,ngas))
        yfn = 0.D0
      END IF
!
      CALL data_exchange(rgpgc)
!
      one  = cte(1.D0)
      DO ijk = 1, ncint
       IF( fl_l(ijk) == 1 ) THEN
         CALL meshinds(ijk,imesh,i,j,k)
         CALL subscr(ijk)
         
         DO ig=1,ngas
	   CALL nb(field,rgpgc(:,ig),ijk)
	   CALL rnb(u,ug,ijk)
	   CALL rnb(w,wg,ijk)

           IF (job_type == '2D') THEN

  	     CALL fsc(yfe(ijk,ig), yft(ijk,ig),      &
                      yfe(imjk,ig), yft(ijkm,ig),   &
                      one, field, u, w, ijk)

           ELSE IF (job_type == '3D') THEN

   	     CALL rnb(v,vg,ijk)
  	     CALL fsc(yfe(ijk,ig), yfn(ijk,ig), yft(ijk,ig),      &
                      yfe(imjk,ig), yfn(ijmk,ig), yft(ijkm,ig),   &
                      one, field, u, v, w, ijk)

           END IF
         END DO
         
       END IF
      END DO
      
      CALL data_exchange(yfe)
      CALL data_exchange(yft)
      IF (job_type == '3D') CALL data_exchange(yfn)

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
           yfx = yfe(ijk,ig) - yfw
	   yfb = yft(ijkm,ig)
	   yfz = yft(ijk,ig) - yfb

           IF (job_type == '3D') THEN
  	     yfs = yfn(ijmk,ig)
  	     yfy = yfn(ijk,ig) - yfs
           END IF

	   rgpgc(ijk,ig) = rgpgcn(ijk,ig)                  &
	                 - dt * indx(i) * yfx              &
	                 - dt * indy(j) * yfy              &
	                 - dt * indz(k) * yfz
 
           IF(rgpgc(ijk,ig) < 0.D0) THEN
             WRITE(8,*) 'Warning!: gas mass is not conserved'
             WRITE(8,128) time, ijk, ig, rgpgc(ijk,ig)
 128         FORMAT('Time= ',F8.3,' Cell= ',I6, ' Specie= ',I2)
             rgpgc(ijk,ig) = 0.D0
           END IF
           
           rgp = rgp + rgpgc(ijk,ig)

         END DO

         ygc(default_gas, ijk) = 1.D0
         DO ig=1,ngas
           IF (ig /= default_gas) THEN
             ygc(ig,ijk)=rgpgc(ijk,ig)/rgp
             ygc(default_gas, ijk) = ygc(default_gas, ijk) - ygc(ig,ijk)
           END IF
         END DO

       END IF
      END DO
!
      DEALLOCATE(yfe, yft)
      IF (job_type == '3D') DEALLOCATE(yfn)
!
      RETURN
      END SUBROUTINE ygas
!----------------------------------------------------------------------
      END MODULE gas_components
!----------------------------------------------------------------------
