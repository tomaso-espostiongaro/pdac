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
      USE convective_fluxes_sc, ONLY: fsc
      USE control_flags, ONLY: job_type
      USE dimensions
      USE domain_decomposition, ONLY: ncint, ncdom, myijk, &
          data_exchange, meshinds
      USE eos_gas, ONLY: rgpgc, rgpgcn, ygc
      USE gas_constants, ONLY: default_gas, gas_type
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE grid, ONLY: dx, dy, dz, indx, indy, indz, inx
      USE grid, ONLY: fl_l
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil, nb, rnb, cte, nb_13
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm
      USE time_parameters, ONLY: dt,time
!
      IMPLICIT NONE
!
      TYPE(stencil) :: one, field, u, v, w
      REAL*8 :: yfw, yfs, yfb, yfx, yfy, yfz
      REAL*8 :: ygcdfg
      REAL*8 :: rgp, rgpinv
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
      CALL data_exchange(rgpgc)
!
      one  = cte(1.D0)

      DO ijk = 1, ncint
       IF( fl_l(ijk) == 1 ) THEN
         CALL meshinds(ijk,imesh,i,j,k)
         CALL subscr(ijk)
         
         DO ig=1,ngas
	   CALL nb_13(field,rgpgc(:,ig),ijk)
           !CALL rnb(u,ug,ijk)
	   !CALL rnb(w,wg,ijk)
           u%w = ug ( imjk )
           u%c = ug ( ijk )
           w%b = wg ( ijkm )
           w%c = wg ( ijk )         
           
           IF (job_type == '2D') THEN

  	     CALL fsc(yfe(ijk,ig), yft(ijk,ig),      &
                      yfe(imjk,ig), yft(ijkm,ig),   &
                      one, field, u, w, ijk)

           ELSE IF (job_type == '3D') THEN
            
	     !CALL rnb(v,vg,ijk)
             v%s = vg ( ijmk )
             v%c = vg ( ijk )
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

	   rgpgc(ijk,ig) = rgpgcn(ijk,ig) - dt * indx(i) * yfx * inx(i)     
           rgpgc(ijk,ig) = rgpgc(ijk,ig) - dt * indy(j) * yfy              
	   rgpgc(ijk,ig) = rgpgc(ijk,ig) - dt * indz(k) * yfz
 
           IF(rgpgc(ijk,ig) < 0.D0) THEN
             WRITE(8,*) 'Warning!: gas mass is not conserved'
             WRITE(8,128) time, ijk, ig 
             WRITE(8,*) rgpgcn(ijk,ig),  rgpgc(ijk,ig)
 128         FORMAT('Time= ',F8.3,' Cell= ',I6, ' Specie= ',I2)
             rgpgc(ijk,ig) = 0.D0
           END IF
           
           rgp = rgp + rgpgc(ijk,ig)

         END DO

         rgpinv = 1.0D0/rgp 

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

       END IF
      END DO
!
      DEALLOCATE(yfe, yft)
      IF (job_type == '3D') DEALLOCATE(yfn)
!
      RETURN
      END SUBROUTINE ygas
!----------------------------------------------------------------------
      SUBROUTINE ygas_3d
!
      USE convective_fluxes_sc, ONLY: fsc, fsc_1st
      USE dimensions
      USE domain_decomposition, ONLY: ncint, ncdom, myijk, &
          data_exchange, meshinds, im1_jp0_kp0_ , ip0_jp0_km1_ , ip0_jm1_kp0_
      USE eos_gas, ONLY: rgpgc, rgpgcn, ygc
      USE flux_limiters, ONLY: muscl
      USE gas_constants, ONLY: default_gas, gas_type
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE grid, ONLY: dx, dy, dz, indx, indy, indz, inx
      USE grid, ONLY: fl_l
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil, nb, rnb, cte, nb_13, first_nb
      USE set_indexes, ONLY: chkstencil
      USE set_indexes, ONLY: subscr, subscr_ygas_1st, subscr_ygas, imjk, ijmk, ijkm
      USE time_parameters, ONLY: dt,time
!
      IMPLICIT NONE
!
      TYPE(stencil) :: one, field, u, v, w
      REAL*8 :: yfw, yfs, yfb, yfx, yfy, yfz
      REAL*8 :: ygcdfg
      REAL*8 :: rgp,rgpinv
      REAL*8 :: rgpgc_tmp
      INTEGER :: i, j, k, imesh
      INTEGER :: ijk
      INTEGER :: ig, dfg, info
      INTEGER :: float_chk
!
      ALLOCATE(yfe(ncdom,ngas), yft(ncdom,ngas))
      yfe = 0.D0; yft = 0.D0
      ALLOCATE(yfn(ncdom,ngas))
      yfn = 0.D0
!
      CALL data_exchange(rgpgc)
!
      one  = cte(1.D0)


      DO ijk = 1, ncint

       IF( fl_l(ijk) == 1 ) THEN

         IF ( muscl == 0) THEN
     
           !CALL meshinds(ijk,imesh,i,j,k)
           CALL subscr_ygas_1st(ijk)
         
           DO ig=1,ngas

             u%w = ug ( imjk )
             u%c = ug ( ijk )
             v%s = vg ( ijmk )
             v%c = vg ( ijk )
             w%b = wg ( ijkm )
             w%c = wg ( ijk )

             !IF( float_chk( u%w ) /= 0 ) WRITE(6,*) 'ygas, wrong u%w: ',u%w
             !IF( float_chk( u%c ) /= 0 ) WRITE(6,*) 'ygas, wrong u%c: ',u%c
             !IF( float_chk( v%s ) /= 0 ) WRITE(6,*) 'ygas, wrong v%s: ',v%s
             !IF( float_chk( v%c ) /= 0 ) WRITE(6,*) 'ygas, wrong v%c: ',v%c
             !IF( float_chk( w%b ) /= 0 ) WRITE(6,*) 'ygas, wrong w%b: ',w%b
             !IF( float_chk( w%c ) /= 0 ) WRITE(6,*) 'ygas, wrong w%c: ',w%c

             CALL first_nb( field, rgpgc(:,ig), ijk )

             !CALL chkstencil( field, info )
             !IF( info /= 0 ) WRITE(6,*) 'ygas, wrong field: ',ijk

             CALL fsc_1st( yfe(ijk,ig), yfn(ijk,ig), yft(ijk,ig),    &
                      yfe(imjk,ig), yfn(ijmk,ig), yft(ijkm,ig),   &
                      one, field, u, v, w)
           
           END DO

         ELSE
	
           !CALL meshinds(ijk,imesh,i,j,k)
           CALL subscr_ygas(ijk)
         
           DO ig=1,ngas

             u%w = ug ( imjk )
             u%c = ug ( ijk )
             v%s = vg ( ijmk )
             v%c = vg ( ijk )
             w%b = wg ( ijkm )
             w%c = wg ( ijk )
  
             CALL nb_13(field,rgpgc(:,ig),ijk)
	     CALL fsc(yfe(ijk,ig), yfn(ijk,ig), yft(ijk,ig),        &
                      yfe(imjk,ig), yfn(ijmk,ig), yft(ijkm,ig),   &
                      one, field, u, v, w, ijk)
           
           END DO

         END IF
         
       END IF
      END DO
      
      CALL data_exchange(yfe)
      CALL data_exchange(yft)
      CALL data_exchange(yfn)

      !    DO ijk = 1, ncint
      !      DO ig = 1, ngas
      !        IF( float_chk( ygc( ig, ijk ) ) /= 0 ) THEN
      !          WRITE(6,*) 'in ygas wrong ygc: ', ijk, ig, ygc(ig,ijk)
      !        END IF
      !        IF( float_chk( yfe( ijk, ig ) ) /= 0 ) THEN
      !          WRITE(6,*) 'in ygas wrong yfe: ', ijk, ig, yfe(ijk,ig)
      !        END IF
      !        IF( float_chk( yft( ijk, ig ) ) /= 0 ) THEN
      !          WRITE(6,*) 'in ygas wrong yft: ', ijk, ig, yft(ijk,ig)
      !        END IF
      !        IF( float_chk( yfn( ijk, ig ) ) /= 0 ) THEN
      !          WRITE(6,*) 'in ygas wrong yfn: ', ijk, ig, yfn(ijk,ig)
      !        END IF
      !      END DO
      !    END DO


      DO ijk = 1, ncint

       IF( fl_l(ijk) == 1 ) THEN

         CALL meshinds(ijk,imesh,i,j,k)

         !CALL subscr(ijk)
         !        Inlined

         imjk   = myijk( im1_jp0_kp0_ , ijk )
         ijmk   = myijk( ip0_jm1_kp0_ , ijk )
         ijkm   = myijk( ip0_jp0_km1_ , ijk )
	   
         yfx = 0.D0
         yfy = 0.D0
         yfz = 0.D0
         rgp = 0.D0

         DO ig = 1, ngas
	   
	   yfw = yfe(imjk,ig)
           yfx = yfe(ijk,ig) - yfw
	   yfb = yft(ijkm,ig)
	   yfz = yft(ijk,ig) - yfb
  	   yfs = yfn(ijmk,ig)
  	   yfy = yfn(ijk,ig) - yfs

	   rgpgc_tmp = rgpgcn(ijk,ig)     
	   rgpgc_tmp = rgpgc_tmp - dt * indx(i) * yfx * inx(i) 
	   rgpgc_tmp = rgpgc_tmp - dt * indy(j) * yfy   
	   rgpgc_tmp = rgpgc_tmp - dt * indz(k) * yfz
           rgpgc(ijk,ig) = rgpgc_tmp
 
           IF( rgpgc(ijk,ig) < 0.D0 ) THEN
             WRITE(8,*) 'Warning!: gas mass is not conserved'
             WRITE(8,128) time, ijk, ig, rgpgc(ijk,ig)
 128         FORMAT('Time= ',F8.3,' Cell= ',I6, ' Specie= ',I2)
             rgpgc(ijk,ig) = 0.D0
           END IF

           IF( float_chk( rgpgc(ijk,ig) ) /= 0 ) &
             WRITE(6,*) 'ygas, wrong rgpgc: ',ijk,ig,rgpgc(ijk,ig)
           
           rgp = rgp + rgpgc(ijk,ig)

         END DO

         ! Carlo Consistency Check 
         IF( rgp >= 1.d-12 ) THEN
           rgpinv = 1.0D0 / rgp
         ELSE
           rgpinv = 0.0d0
         END IF

         IF( float_chk( rgpinv ) /= 0 ) WRITE(6,*) 'ygas, wrong rgpinv: ',ijk,ig,rgpinv,rgp
         
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

       END IF

      END DO
!
      DEALLOCATE(yfe, yft)
      DEALLOCATE(yfn)
!
      RETURN
      END SUBROUTINE ygas_3d
!----------------------------------------------------------------------
      END MODULE gas_components
!----------------------------------------------------------------------
