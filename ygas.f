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
      USE convective_fluxes, ONLY: fsc
      USE eos_gas, ONLY: rgpgc, rgpgcn, ygc
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE grid, ONLY: dx, dy, dz, indx, indy, indz
      USE grid, ONLY: ncint, ncdom, myijk, data_exchange
      USE grid, ONLY: fl_l
      USE set_indexes
      USE time_parameters, ONLY: dt
      USE indijk_module, ONLY: ip0_jp0_kp0_
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
      DO ijk = 1, ncint
       imesh = myijk( ip0_jp0_kp0_, ijk)
       IF( fl_l(ijk) == 1 ) THEN
         i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
         j = MOD( imesh - 1, nx*ny ) / nx + 1
         k = ( imesh - 1 ) / ( nx*ny ) + 1
         CALL subscr(ijk)
         
         DO ig=1,ngas
	   one  = cte(1.D0)
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

      rgp = 0.D0
      DO ijk = 1, ncint
       imesh = myijk( ip0_jp0_kp0_, ijk)
       IF( fl_l(ijk) == 1 ) THEN
         i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
         j = MOD( imesh - 1, nx*ny ) / nx + 1
         k = ( imesh - 1 ) / ( nx*ny ) + 1
         CALL subscr(ijk)
	   
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

           IF(rgpgc(ijk,ig) < 0.D0) CALL error('ygas','gas mass not conserved',1)
           rgp = rgp + rgpgc(ijk,ig)

         END DO
       END IF
      END DO
	
      DO ig=1,ngas
        ygc(ig,ijk)=rgpgc(ijk,ig)/rgp
      END DO
!
      DEALLOCATE(yfe, yfn, yft)
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE gas_components
!----------------------------------------------------------------------
