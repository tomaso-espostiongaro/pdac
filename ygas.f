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
      USE dimensions
      USE convective_fluxes_sc, ONLY: fsc
      USE eos_gas, ONLY: rgpgc, rgpgcn, ygc
      USE gas_constants, ONLY: default_gas
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE grid, ONLY: dx, dy, dz, indx, indy, indz
      USE grid, ONLY: ncint, ncdom, myijk, data_exchange
      USE grid, ONLY: fl_l
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
      ALLOCATE(yfe(ncdom,ngas), yfn(ncdom,ngas), yft(ncdom,ngas))
      yfe = 0.D0; yfn = 0.D0; yft = 0.D0
!
      CALL data_exchange(rgpgc)
!
      one  = cte(1.D0)
      DO ijk = 1, ncint
       IF( fl_l(ijk) == 1 ) THEN
         imesh = myijk( ip0_jp0_kp0_, ijk)
         i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
         j = MOD( imesh - 1, nx*ny ) / nx + 1
         k = ( imesh - 1 ) / ( nx*ny ) + 1
         CALL subscr(ijk)
         
         DO ig=1,ngas
	   CALL nb(field,rgpgc(:,ig),ijk)
	   CALL rnb(u,ug,ijk)
	   CALL rnb(v,vg,ijk)
	   CALL rnb(w,wg,ijk)

	   CALL fsc(yfe(ijk,ig), yfn(ijk,ig), yft(ijk,ig),      &
                    yfe(imjk,ig), yfn(ijmk,ig), yft(ijkm,ig),   &
                    one, field, u, v, w, ijk)
         END DO
         
       END IF
      END DO
      
      CALL data_exchange(yfe)
      CALL data_exchange(yfn)
      CALL data_exchange(yft)

      DO ijk = 1, ncint
       IF( fl_l(ijk) == 1 ) THEN
         imesh = myijk( ip0_jp0_kp0_, ijk)
         i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
         j = MOD( imesh - 1, nx*ny ) / nx + 1
         k = ( imesh - 1 ) / ( nx*ny ) + 1
         CALL subscr(ijk)
	   
         rgp = 0.D0
         DO ig=1,ngas
	   
	   yfw = yfe(imjk,ig)
	   yfs = yfn(ijmk,ig)
	   yfb = yft(ijkm,ig)

           yfx = yfe(ijk,ig) - yfw
	   yfy = yfn(ijk,ig) - yfs
	   yfz = yft(ijk,ig) - yfb

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
      DEALLOCATE(yfe, yfn, yft)
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE gas_components
!----------------------------------------------------------------------
