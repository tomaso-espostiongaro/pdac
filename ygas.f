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
      ALLOCATE(yfe(ngas,ncdom), yfn(ngas,ncdom), yft(ngas,ncdom))
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
	   field = nb(rgpgc,ig,ijk)
	   u    = rnb(ug,ijk)
	   v    = rnb(vg,ijk)
	   w    = rnb(wg,ijk)
	   CALL fsc(yfe(ig,ijk), yfn(ig,ijk), yft(ig,ijk),      &
                    yfe(ig,imjk), yfn(ig,ijmk), yft(ig,ijkm),   &
                    one, field, u, v, w, ijk)
         END DO
         
       END IF
      END DO
      
      CALL data_exchange(yfe)
      CALL data_exchange(yfn)
      CALL data_exchange(yft)

      DO ijk = 1, ncint
       imesh = myijk( ip0_jp0_kp0_, ijk)
       IF( fl_l(ijk) == 1 ) THEN
         i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
         j = MOD( imesh - 1, nx*ny ) / nx + 1
         k = ( imesh - 1 ) / ( nx*ny ) + 1
         CALL subscr(ijk)
	   
         DO ig=1,ngas
	   
	   yfw = yfe(ig,imjk)
	   yfs = yfn(ig,ijmk)
	   yfb = yft(ig,ijkm)

           yfx = yfe(ig,ijk) - yfw
	   yfy = yfn(ig,ijk) - yfs
	   yfz = yft(ig,ijk) - yfb

	   rgpgc(ig,ijk) = rgpgcn(ig,ijk)                  &
	                 - dt * indx(i) * yfx              &
	                 - dt * indy(j) * yfy              &
	                 - dt * indz(k) * yfz

           IF(rgpgc(ig,ijk) < 0.D0) CALL error('ygas','gas mass not conserved',1)

         END DO
       END IF
      END DO
	
      rgp = SUM(rgpgc(:,ijk))
      ygc(:,ijk)=rgpgc(:,ijk)/rgp
!
      DEALLOCATE(yfe, yfn, yft)
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE gas_components
!----------------------------------------------------------------------
