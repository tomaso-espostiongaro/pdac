!-----------------------------------------------------------------------
      MODULE boundary_conditions
!-----------------------------------------------------------------------
!
      USE atmospheric_conditions, ONLY: gravx, gravy, gravz
      USE dimensions
      USE gas_solid_density, ONLY: rgp, rlk, rog
      USE gas_solid_temperature, ONLY: sieg, tg, sies, ts
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE grid, ONLY: dx, dy, dz, xb, zb
      USE eos_gas, ONLY: xgc, ygc
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p, ep
      USE time_parameters, ONLY: dt, time, sweep

      IMPLICIT NONE

      INTEGER :: n0, n1, n2

      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE boundary
!
! ... This routine computes (x,y,z) boundary conditions 
!
      USE control_flags, ONLY: job_type, lpr
      USE io_files, ONLY: tempunit
      USE domain_decomposition, ONLY: ncint, myijk, meshinds
      USE grid, ONLY: flag, x, y, z, xb, yb, zb
      USE grid, ONLY: slip_wall, noslip_wall, fluid, int_immb, ext_immb
      USE grid, ONLY: free_io, nrfree_io, inlet_cell, vent_cell
      USE immersed_boundaries, ONLY: fptx, fpty, fptz
      USE immersed_boundaries, ONLY: numx, numy, numz, immb
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE parallel, ONLY: mpime, root
      USE set_indexes, ONLY: subscr
      USE set_indexes, ONLY: ipjk, imjk, ippjk, immjk, ijpk, ipjpk,    &
        imjpk, ijmk, ipjmk, imjmk, ijppk, ijmmk, ijkp, ipjkp, imjkp,   &
        ijpkp, ijmkp, ijkm, ipjkm, imjkm, ijpkm, ijmkm, ijkpp, ijkmm
      USE vent_conditions, ONLY: update_ventc, random_switch, irand, &
                                 update_vent_cell
!
      IMPLICIT NONE
!
      INTEGER :: ijk, i, j, k, imesh, ig, is
      INTEGER :: fp, np
      REAL*8 :: d1, d2 
      REAL*8 :: vel
      INTEGER :: fx, fy, fz
      LOGICAL :: forced
!
      IF (irand >= 1) CALL random_switch(sweep)
!
      mesh_loop: DO ijk = 1, ncint

        fx = 0
        fy = 0 
        fz = 0
        forced = .FALSE.

        CALL subscr(ijk)
        CALL meshinds(ijk,imesh,i,j,k)

!
! ... Update inlet cells for non-stationnary boundary conditions
!
        IF (irand >= 1) THEN
          IF (flag(ijk) == inlet_cell) THEN
            CALL update_vent_cell(ijk)
          ELSE IF (flag(ijk) == vent_cell) THEN
            CALL update_ventc(ijk,imesh,sweep)
          END IF
        END IF
!
! ... In fluid cells and immersed boundaries, update the 
! ... neighbours on boundaries
!
        IF( flag(ijk) == int_immb .OR. flag(ijk) == ext_immb) THEN
          
          ! ... Check if 'ijk' is a forcing point
          !
          fx = numx(ijk)
          IF (job_type == '3D') fy = numy(ijk)
          fz = numz(ijk)
          forced = (fx/=0 .OR. fy/=0 .OR. fz/=0)
          IF (.NOT.forced) CALL error('boundary','control forcing points',1)
          
          ! ... Compute the pseudo-velocities
          ! ... that are used in the "immersed boundary" technique ...
          !
          IF (job_type == '2D') THEN

            IF( fx/=0 ) THEN
              vel = velint(fptx(fx), ug, ijk, xb, y, z)  
              fptx(fx)%vel = vel
              !
              ! ... Set the pressure in non-resolved forcing points
              ! ... (zero-gradient)
              IF( fptx(fx)%int == 17 ) THEN
                      IF (flag(ijk) == fluid) THEN
                              p(ipjk) = p(ijk)
                      ELSE
                              p(ijk) = p(ipjk)
                      END IF
              END IF
              !
              ! ... Initialize x-velocity in the forced points
              ug(ijk) = vel
              us(ijk,:) = vel
            END IF

            IF( fz/=0 ) THEN
              vel = velint(fptz(fz), wg, ijk, x, y, zb)  
              fptz(fz)%vel = vel
              !
              ! ... Set the pressure in non-resolved forcing points
              IF( fptz(fz)%int == 17 ) THEN
                      p(ijk) = p(ijkp)
              END IF
              !
              ! ... Initialize z-velocity in the forced points
              wg(ijk) = vel
              ws(ijk,:) = vel
            END IF

          ELSE IF (job_type == '3D') THEN

            IF( fx/=0 ) THEN
              vel = velint3d(fptx(fx), ug, ijk, xb, y, z)  
              fptx(fx)%vel = vel
              !
              ! ... Set the pressure in non-resolved forcing points
              IF( fptx(fx)%int == 17 ) THEN
                      IF (flag(ijk) == fluid) THEN
                              p(ipjk) = p(ijk)
                      ELSE
                              p(ijk) = p(ipjk)
                      END IF
              END IF
              !
              ! ... Initialize x-velocity in the forced points
              ug(ijk) = vel
              us(ijk,:) = vel
            END IF
            
            IF( fy/=0 ) THEN
              vel = velint3d(fpty(fy), vg, ijk, x, yb, z)  
              fpty(fy)%vel = vel

              ! ... Set the pressure in non-resolved forcing points
              IF( fpty(fy)%int == 17 ) THEN
                      IF (flag(ijk) == fluid) THEN
                              p(ijpk) = p(ijk)
                      ELSE
                              p(ijk) = p(ijpk)
                      END IF
              END IF
              
              ! ... Initialize y-velocity in the forced points
              vg(ijk) = vel
              vs(ijk,:) = vel
            END IF
            
            IF( fz/=0 ) THEN
              vel = velint3d(fptz(fz), wg, ijk, x, y, zb)  
              fptz(fz)%vel = vel

              ! ... Set the pressure in non-resolved forcing points
              IF( fptz(fz)%int == 17 ) THEN
                      p(ijk) = p(ijkp)
              END IF
              
              ! ... Initialize z-velocity in the forced points
              wg(ijk) = vel
              ws(ijk,:) = vel
            END IF
          END IF
!
! ... In fluid cells update the neighbours on boundaries
!
        ELSE IF( flag(ijk) == fluid ) THEN
!
! ***** East boundary conditions ***** !
!
            n2   = ipjk
	    n1   = ijk
	    n0   = imjk

            SELECT CASE ( flag( n2 ) ) 
!
            CASE (slip_wall) 
!
              wg(n2)   = wg(n1)
              ws(n2,:) = ws(n1,:)

	      IF (job_type == '3D') THEN
                vg(n2)   = vg(n1)
                vs(n2,:) = vs(n1,:)
              END IF
!
            CASE (noslip_wall)
!
              ug(n2)   = 0.D0
              us(n2,:) = 0.D0
              
              IF(flag(ipjkp) /= fluid) THEN
                wg(n2)   = -wg(n1)
                DO is = 1, nsolid
                  IF (rlk(ijk,is) > 0.D0) ws(n2,is) = -ws(n1,is)
                END DO
              END IF

              IF (flag(ipjpk) /= fluid .AND. job_type == '3D') THEN
                vg(n2)   = -vg(n1)
                DO is = 1, nsolid
                  IF (rlk(ijk,is) > 0.D0) vs(n2,:) = -vs(n1,:)
                END DO
              END IF
!
            CASE (free_io)
!
              d1 = dx(i)
	      d2 = dx(i+1)

              ! ... Compute the normal component of the velocities 
              ! ... and scalar fields 
              !
              CALL ent_inout4( ug(n0), ug(n1), ug(n2),              &
                      us(n0,:), us(n1,:), us(n2,:), d1, d2, gravx, k )

              ! ... Compute tangential components of velocities
              !
              IF(flag(ipjkp) == free_io) THEN
                wg(n2) = wg(n1)
        	ws(n2,:)=ws(n1,:)
              END IF
	      IF (job_type == '3D' .AND. (flag(ipjpk) == free_io)) THEN
                vg(n2) = vg(n1)
        	vs(n2,:)=vs(n1,:)
	      END IF
!
            CASE (nrfree_io)
!
              d1 = dx(i)
	      d2 = dx(i+1)

              ! ... Compute the normal component of the velocities 
              ! ... and scalar fields
              !
              CALL ent_inout6( ug(n1), ug(n2), us(n1,:), us(n2,:), d1, d2, k)

              ! ... Compute tangential components of velocities
              !
              IF(flag(ipjkp) == nrfree_io) THEN
                wg(n2) = wg(n1)
        	ws(n2,:)=ws(n1,:)
              END IF
	      IF (job_type == '3D' .AND. (flag(ipjpk) == nrfree_io)) THEN
                vg(n2) = vg(n1)
        	vs(n2,:)=vs(n1,:)
	      END IF
!              
            CASE DEFAULT
!
              CONTINUE

            END SELECT
!
! ***** West boundary conditions ***** !
!
            n2 = imjk
	    n1 = ijk
	    n0 = ipjk

            SELECT CASE ( flag( n2 ) )
!
            CASE (slip_wall)
!
              wg(n2)   = wg(n1)
              ws(n2,:) = ws(n1,:)
	      IF (job_type == '3D') THEN
                vg(n2)   = vg(n1)
                vs(n2,:) = vs(n1,:)
	      END IF
!
            CASE (noslip_wall)

              ug(n2)   = 0.D0
              us(n2,:) = 0.D0

              IF ( flag(imjkp) /= fluid ) THEN
                wg(n2)   = -wg(n1)
                DO is = 1, nsolid
                  IF (rlk(ijk,is) > 0.D0) ws(n2,is) = -ws(n1,is)
                END DO
	      END IF

              IF ( flag(imjkp) /= fluid .AND. job_type == '3D') THEN
                vg(n2)   = -vg(n1)
                DO is = 1, nsolid
                  IF (rlk(ijk,is) > 0.D0) vs(n2,is) = -vs(n1,is)
                END DO
	      END IF
!
            CASE (free_io)
!	    
              d1 = dx(i)
              d2 = dx(i-1)

              ! ... Compute the normal component of the velocities 
              ! ... and scalar fields
              !
              CALL wsb_inout4( ug(n1), ug(n2), us(n1,:), us(n2,:), d1, d2, k )

              ! ... Compute tangential components of velocities
              !
              IF(flag(imjkp) == free_io) THEN
                wg(n2)   = wg(n1)
        	ws(n2,:) = ws(n1,:)
              END IF
	      IF (job_type == '3D' .AND. (flag(imjpk) == free_io)) THEN
                vg(n2)   = vg(n1)
        	vs(n2,:) = vs(n1,:)
	      END IF
!              
            CASE (nrfree_io)
!	    
              d1 = dx(i)
              d2 = dx(i-1)

              ! ... Compute the normal component of the velocities 
              ! ... and scalar fields
              !
              CALL wsb_inout6( ug(n1), ug(n2), us(n1,:), us(n2,:), d1, d2, k )

              ! ... Compute tangential components of velocities
              !
              IF(flag(imjkp) == nrfree_io) THEN
                wg(n2)   = wg(n1)
        	ws(n2,:) = ws(n1,:)
              END IF
	      IF (job_type == '3D' .AND. (flag(imjpk) == nrfree_io)) THEN
                vg(n2)   = vg(n1)
        	vs(n2,:) = vs(n1,:)
	      END IF
!              
            CASE DEFAULT
!
              CONTINUE

            END SELECT

            IF (job_type == '3D') THEN
!
! ***** North boundary conditions ***** !
!
              n2   = ijpk
  	      n1   = ijk
  	      n0   = ijmk
  
              SELECT CASE (  flag( n2 ) )
!  
              CASE (slip_wall)
!  
                ug(n2)   = ug(n1)
                us(n2,:) = us(n1,:)
                wg(n2)   = wg(n1)
                ws(n2,:) = ws(n1,:)
!  
              CASE (noslip_wall)
!  
                vg(n2)   = 0.D0
                vs(n2,:) = 0.D0
                
                IF(flag(ipjpk) /= fluid) THEN
                  ug(n2)   = -ug(n1)
                  DO is = 1, nsolid
                    IF (rlk(ijk,is) > 0.D0) us(n2,is) = -us(n1,is)
                  END DO
                END IF

                IF(flag(ijpkp) /= fluid) THEN
                  wg(n2)   = -wg(n1)
                  DO is = 1, nsolid
                    IF (rlk(ijk,is) > 0.D0) ws(n2,is) = -ws(n1,is)
                  END DO
                END IF
!  
              CASE (free_io)
!  
                d1 = dy(j)
  	        d2 = dy(j+1)
          
                ! ... Compute the normal component of the velocities 
                ! ... and scalar fields
                !
                CALL ent_inout4( vg(n0), vg(n1), vg(n2),              &
                             vs(n0,:), vs( n1,:), vs(n2,:), d1, d2, gravy, k )

                ! ... Compute tangential components of velocities
                !
                IF(flag(ijpkp) == free_io) THEN
                  wg(n2)   = wg(n1)
        	  ws(n2,:) = ws(n1,:)
                END IF
	        IF(flag(ipjpk) == free_io) THEN
                  ug(n2)   = ug(n1)
          	  us(n2,:) = us(n1,:)
	        END IF
!
              CASE (nrfree_io)
!  
                d1 = dy(j)
  	        d2 = dy(j+1)
          
                ! ... Compute the normal component of the velocities 
                ! ... and scalar fields
                !
                CALL ent_inout6( vg(n1), vg(n2), vs(n1,:), vs(n2,:), d1, d2, k)

                ! ... Compute tangential components of velocities
                !
                IF(flag(ijpkp) == nrfree_io) THEN
                  wg(n2)   = wg(n1)
        	  ws(n2,:) = ws(n1,:)
                END IF
	        IF(flag(ipjpk) == nrfree_io) THEN
                  ug(n2)   = ug(n1)
          	  us(n2,:) = us(n1,:)
	        END IF
!
              CASE DEFAULT
!
                CONTINUE

              END SELECT

!
! ***** South boundary conditions ***** !
!
              n2 = ijmk
              n1 = ijk
              n0 = ijpk

              SELECT CASE (  flag( n2 ) )
!
              CASE (slip_wall)
!
                ug(n2)   = ug(n1)
                us(n2,:) = us(n1,:)
                wg(n2)   = wg(n1)
                ws(n2,:) = ws(n1,:)
!  
              CASE (noslip_wall)
!
                vg(n2)   = 0.D0
                vs(n2,:) = 0.D0

                IF(flag(ipjmk) /= fluid) THEN
                  ug( n2 ) = -ug( n1 )
                  DO is = 1, nsolid
                    IF (rlk(ijk,is) > 0.D0) us(n2,is) = -us(n1,is)
                  END DO
                  
                END IF
                IF(flag(ijmkp) /= fluid) THEN
                  wg( n2 ) = -wg( n1 )
                  DO is = 1, nsolid
                    IF (rlk(ijk,is) > 0.D0) ws(n2,is) = -ws(n1,is)
                  END DO
                END IF
!
              CASE (free_io)
!  	    
                d1 = dy(j)
                d2 = dy(j-1)

                ! ... Compute the normal component of the velocities 
                ! ... and scalar fields
                !
                CALL wsb_inout4( vg(n1), vg(n2), vs(n1,:), vs(n2,:), d1, d2, k )

                ! ... Compute tangential components of velocities
                !             
                IF(flag(ijmkp) == free_io) THEN
                  wg(n2)   = wg(n1)
        	  ws(n2,:) = ws(n1,:)
                END IF
	        IF(flag(ipjmk) == free_io) THEN
                  ug(n2)   = ug(n1)
        	  us(n2,:) = us(n1,:)
	        END IF
!
              CASE (nrfree_io)
!  	    
                d1 = dy(j)
                d2 = dy(j-1)

                ! ... Compute the normal component of the velocities 
                ! ... and scalar fields
                !
                CALL wsb_inout6( vg(n1), vg(n2), vs(n1,:), vs(n2,:), d1, d2, k )

                ! ... Compute tangential components of velocities
                !             
                IF(flag(ijmkp) == nrfree_io) THEN
                  wg(n2)   = wg(n1)
        	  ws(n2,:) = ws(n1,:)
                END IF
	        IF(flag(ipjmk) == nrfree_io) THEN
                  ug(n2)   = ug(n1)
        	  us(n2,:) = us(n1,:)
	        END IF
!
              CASE DEFAULT
!
                CONTINUE

              END SELECT

            END IF

!
! ***** Top boundary conditions ***** !
!
            n2   =  ijkp
            n1   =  ijk
            n0   =  ijkm

            SELECT CASE ( flag( n2 ) )
!
            CASE (slip_wall)
!
              ug(n2)   = ug(n1)
              us(n2,:) = us(n1,:)
	      IF (job_type == '3D') THEN
                vg(n2)   = vg(n1)
                vs(n2,:) = vs(n1,:)
	      END IF
!
            CASE (noslip_wall)
!  
              wg(n2)   = 0.D0
              ws(n2,:) = 0.D0

              IF (flag(ipjkp) /= fluid) THEN
                ug(n2)   = -ug(n1)
                DO is = 1, nsolid
                  IF (rlk(ijk,is) > 0.D0) us(n2,is) = -us(n1,is)
                END DO
              END IF

	      IF (flag(ijpkp) /= fluid .AND. job_type == '3D') THEN
                vg(n2)   = -vg(n1)
                DO is = 1, nsolid
                  IF (rlk(ijk,is) > 0.D0) vs(n2,is) = -vs(n1,is)
                END DO
	      END IF
!
            CASE (free_io)
!
              d1 = dz(k)
	      d2 = dz(k+1)

              ! ... Compute the normal component of the velocities 
              ! ... and scalar fields
              !
              CALL ent_inout4( wg(n0), wg(n1), wg(n2),              &
                             ws(n0,:), ws(n1,:), ws(n2,:), d1, d2, gravz, k )

              ! ... Compute tangential components of velocities
              !
              IF(flag(ipjkp) == free_io) THEN
                 ug(n2)   = ug(n1)
                 us(n2,:) = us(n1,:)
              END IF
              IF(flag(ijpkp) == free_io .AND. job_type == '3D') THEN
                 vg(n2)   = vg(n1)
                 vs(n2,:) = vs(n1,:)
              END IF
!
            CASE (nrfree_io)
!
              d1 = dz(k)
              d2 = dz(k+1)

              ! ... Compute the normal component of the velocities 
              ! ... and scalar fields
              !
              CALL ent_inout6( wg(n1), wg(n2), ws(n1,:), ws(n2,:), d1, d2, k)
                             
              ! ... Compute tangential components of velocities
              !
              IF(flag(ipjkp) == nrfree_io) THEN
                 ug(n2)   = ug(n1)
                 us(n2,:) = us(n1,:)
              END IF
              IF(flag(ijpkp) == nrfree_io .AND. job_type == '3D') THEN
                 vg(n2)   = vg(n1)
                 vs(n2,:) = vs(n1,:)
              END IF
!
            CASE DEFAULT
!
              CONTINUE

            END SELECT

            ! ... set upper corners velocities
            !
            IF ( k == (nz-1) ) THEN
              IF( ( i == (nx-1) ) .AND. ( j == (ny-1) ) ) THEN
                ! ug(ipjp) = ug(ipj)
              ELSE IF( ( i == 2 ) .AND. ( j == (ny-1) ) ) THEN
              ELSE IF( ( i == (nx-1) ) .AND. ( j == 2 ) ) THEN
                ! ug(imjp) = ug(imj)
              ENDIF
            END IF
!
! ***** Bottom boundary conditions ***** !
!
            n1 = ijk
            n2 = ijkm

            SELECT CASE ( flag( n2 ) )
!
            CASE (slip_wall)
!
              IF(flag(ipjkm) == slip_wall) THEN 
                ug( n2 ) = ug( n1 )
                us(n2,:) = us(n1,:)
              END IF

              IF (flag(ijpkm) == slip_wall .AND. job_type == '3D') THEN
                vg( n2 ) = vg( n1 )
                vs(n2,:) = vs(n1,:)
              END IF
!
            CASE (noslip_wall)
!
              wg(n2)   = 0.D0
              ws(n2,:) = 0.D0
                
              IF(flag(ipjkm) == noslip_wall) THEN 
                ug(n2)   = -ug(n1)
                DO is = 1, nsolid
                  IF (rlk(ijk,is) > 0.D0) us(n2,is) = -us(n1,is)
                END DO
              END IF

              IF(flag(ijpkm) == noslip_wall .AND. job_type == '3D') THEN
                vg(n2)   = -vg(n1)
                DO is = 1, nsolid
                  IF (rlk(ijk,is) > 0.D0) vs(n2,is) = -vs(n1,is)
                END DO
              END IF
!
            CASE DEFAULT
!
              CONTINUE

            END SELECT

            ! ... Set lower corners velocities
            !
!            IF (k == 2) THEN
!              IF( ( i == (nx-1) ) .AND. ( j == (ny-1) ) ) THEN
!                ! ug(ipjm) = ug(ipj)
!              ELSE IF( ( i == 2 ) .AND. ( j == (ny-1) ) ) THEN
!              ELSE IF( ( i == (nx-1) ) .AND. ( j == 2 ) ) THEN
!              ELSE IF( ( i == 2 ) .AND. ( j == 2 ) ) THEN
!                ! ug(imjm) = ug(imj)
!              ENDIF
!            END IF

        END IF
      END DO mesh_loop
!
      IF (lpr > 1 .AND. immb >= 1) THEN
        IF (mpime == root) THEN
          OPEN(UNIT=tempunit,FILE='fptx.dat',STATUS='UNKNOWN')
          DO np = 1, SIZE(fptx)
            WRITE(tempunit,33) np, fptx(np)
          END DO
          CLOSE(tempunit)
          IF (job_type == '3D') THEN
            OPEN(UNIT=tempunit,FILE='fpty.dat',STATUS='UNKNOWN')
            DO np = 1, SIZE(fpty)
              WRITE(tempunit,33) np, fpty(np)
            END DO
            CLOSE(tempunit)
          END IF
          OPEN(UNIT=tempunit,FILE='fptz.dat',STATUS='UNKNOWN')
          DO np = 1, SIZE(fptz)
            WRITE(tempunit,33) np, fptz(np)
          END DO
          CLOSE(tempunit)
        END IF
 33   FORMAT(5(I6),10(F16.3))
      END IF
!
      RETURN
      END SUBROUTINE boundary
!----------------------------------------------------------------------
!     I N - O U T   F L O W   P R O C E D U R E S 
!----------------------------------------------------------------------
      SUBROUTINE ent_inout4(umn, ucn, upn, usmn, uscn, uspn, d1, d2, grav, k)
!
! ... This routine computes the free in/outflow conditions in the boundary
! ... cell, i.e. the normal component of the velocity and the scalar fields
! ... for /// E a s t , N o r t h , T o p /// boundaries

      USE atmospheric_conditions, ONLY: p_atm, t_atm
      USE gas_constants, ONLY: gmw, gammaair, gamn, rgas
      USE gas_constants, ONLY: gas_type

      REAL*8, INTENT(IN) :: umn, ucn
      REAL*8, INTENT(INOUT) :: upn
      REAL*8, INTENT(IN) :: usmn(:), uscn(:)
      REAL*8, INTENT(INOUT) :: uspn(:)
      REAL*8, INTENT(IN) :: d1, d2, grav
      INTEGER, INTENT(IN) :: k

      REAL*8 :: prif, pnn2, p1nn
      REAL*8 :: zrif, trif, rhorif, cost, costc, rhorifm1
      REAL*8 :: rmcn, rmcnn, rmmn, rm1nn, rm1knn, rm2n, rm1n, rm0n
      REAL*8 :: u1n, u2n, ucnn, upnn
      REAL*8 :: eps, epc, epcn, ep1nn, epnn
      REAL*8 :: tcn, tmn, t1nn, t1n, t0n, t2n
      REAL*8 :: mg
      REAL*8 :: dc
      REAL*8 :: ind1,ind2,indc
      REAL*8 :: pf1, pf2, pfd

      INTEGER :: ig, is
!
! ... Definitions
!
      u1n = ( umn + ucn ) * 0.5D0
      u2n = ( upn + ucn ) * 0.5D0
      dc = ( d2  + d1 ) * 0.5D0
      ind1 = 1.0D0 / d1
      ind2 = 1.0D0 / d2
      indc = 1.0D0 / dc
!
! ... interpolation of mid-point values

      t0n  = tg( n0 )
      t1n  = tg( n1 )
      t2n  = tg( n2 )
      tmn  = ( t0n + t1n ) * 0.5D0
      tcn  = ( t1n + t2n ) * 0.5D0
      epcn = ( ep( n1 ) + ep( n2 ) ) * 0.5D0

      IF( ucn > 0.D0 ) THEN

! ...  OUTFLOW ...
!
! ...  Extrapolations
       ucnn = u1n
       upn = u2n
       !
       t1nn = tmn
!
! ... calculation of the gas volumetric fraction at time (n+1)dt
! ... from the mass balance equation of solids in cell (ijk)

        epc = 0.D0
        DO is = 1, nsolid
          eps = rlk(n1,is) * inrl(is) - ( dt * inrl(is) * ind1 )  *          &
                ( uscn(is)*rlk(n1,is) - usmn(is)*rlk(n0,is) )
          epc = epc + eps
        END DO
        ep1nn = 1.D0 - epc
!
! ... calculation of the mixture density 
! ... from a mass balance equation for the mixture.
! ... Use gas velocity for mixture.

        rm2n=0.D0
        rm1n=0.D0
        rm0n=0.D0
        DO is=1,nsolid
          rm2n=rm2n+rlk(n2,is)
          rm1n=rm1n+rlk(n1,is)
          rm0n=rm0n+rlk(n0,is)
        END DO
        rm1knn = rm1n
        rm2n = rm2n + rgp(n2)
        rm1n = rm1n + rgp(n1)
        rm0n = rm0n + rgp(n0)
        rmcn=(rm2n+rm1n)*0.5D0
        rmmn=(rm1n+rm0n)*0.5D0
!
        rmcnn = rmcn-dt*indc*(u2n*rmcn-u1n*rmmn)
!
! ... Calculation of fluid pressure from the gas equation of state and 
! ... the mixture density transport equation
!
        mg = 0.D0
        DO ig = 1, ngas
          mg = mg + xgc(n1,ig) * gmw( gas_type( ig ) )
        END DO
        rm1nn = ( rgas * t1nn ) / ( ep1nn * mg )      
        p1nn  = rm1nn * (- rm1knn + rm1n - dt * ind1 * ( rm1n * ucn - rm0n * umn ) )

        IF (epc < 1.0D-8) p1nn = p(n1)
!
! ... Calculation of the advanced-time fluid pressure from Momentum balance
! ... equation of the mixture

        p(n2) = -rmcnn*ucnn + rmcn*ucn - dt*indc*( u2n*ucn*rmcn - u1n*umn*rmmn) 
        p(n2) = p(n2) + grav * dt * rmcn
        p(n2) = p(n2) / ( dt * indc ) + p1nn
!
        ep(n2) = ep(n1)
        tg(n2) = tg(n1)
        DO ig = 1, ngas
           ygc(n2,ig) = ygc(n2,ig) - &
                          ucn * dt * ind2 * ( ygc(n2,ig) - ygc(n1,ig) )
        END DO
!
! ... Correct non-physical pressure
!
        IF ( p(n2) <= 0.0D0 ) p(n2) = p(n1)
!
      ELSE IF ( ucn <  0.D0 ) THEN

! ... INFLOW ...
!
! ... Extrapolation
        upn = ucn
!
        zrif = zb(k) + 0.5D0 * ( dz(1) - dz(k) )  ! DOMANDA perche dz(1)
        prif = p_atm(k)
        trif = t_atm(k)

        rhorifm1 = rgas * trif / ( prif * gmw(6) )   !  rhorif ** (-1)
        cost  = prif ** ( 1.D0 / gammaair ) * rhorifm1
        costc = ( gammaair - 1.D0 ) / ( gammaair * cost )
!
! ... Adiabatic inflow
 
        pf1 = prif**gamn
        pf2 = u2n * u2n * 0.5D0 * costc
        pfd = pf1 - pf2

        IF( pfd >= 0 ) THEN
          p(n2) = pfd ** ( 1.d0 / gamn )
        ELSE
          p(n2) = prif
        END IF
!
        ep(n2) = 1.D0
        tg(n2) = trif
        DO ig = 1, ngas
          ygc(n2,ig) = ygc(n1,ig)
        END DO

! ... Correct non-physical pressure

        IF ( p(n2) < 0.0D0 ) p(n2) = prif

!
      ELSE IF ( ucn == 0.D0 ) THEN
!
        upn      = ucn
        p(n2)    = p(n1)
        ep(n2)   = ep(n1)
        tg(n2)   = tg(n1)
        DO ig=1,ngas
          ygc(n2,ig) = ygc(n1,ig)
        END DO

      ENDIF
!
! ... Set primary variables
!
      DO is = 1, nsolid
        IF( uscn(is) >= 0.D0 ) THEN
          rlk(n2,is) = rlk(n1,is)
          uspn(is)  = uscn(is)
        ELSE
          rlk(n2,is) = 0.0D0
          uspn(is)  = 0.0D0
        ENDIF
      END DO
 
      sieg(n2) = sieg(n1)
      DO is = 1, nsolid
        sies(n2,is) = sies(n1,is)
      END DO
!                
      RETURN
      END SUBROUTINE ent_inout4
!----------------------------------------------------------------------
      SUBROUTINE wsb_inout4(ucn, umn, uscn, usmn, d1, d2, k)
!
! ... This routine computes the free in/outflow conditions in the boundary
! ... cell, i.e. the normal component of the velocity and the scalar fields
! ... for /// W e s t,  S o u t h ///  boundaries

      USE atmospheric_conditions, ONLY: p_atm, t_atm
      USE gas_constants, ONLY: gmw, gammaair, gamn, rgas
      USE gas_constants, ONLY: gas_type

      REAL*8, INTENT(IN) :: ucn
      REAL*8, INTENT(INOUT) :: umn
      REAL*8, INTENT(IN) :: uscn(:)
      REAL*8, INTENT(INOUT) :: usmn(:)
      REAL*8, INTENT(IN) :: d1, d2
      INTEGER, INTENT(IN) :: k

      REAL*8 :: prif, pnn2, p1nn
      REAL*8 :: zrif, trif, rhorif, cost, costc
      REAL*8 :: rmcn, rmcnn, rmpn, rm1nn, rm1knn, rm2n, rm1n, rm0n
      REAL*8 :: u1n, u2n, umnn
      REAL*8 :: eps, epc, epcn, ep1nn
      REAL*8 :: tcn, tpn, t1nn, t1n, t0n, t2n
      REAL*8 :: mg
      REAL*8 :: dc
      REAL*8 :: ind1,ind2,indc
      REAL*8 :: pf1, pf2, pfd

      INTEGER :: ig, is
!
! ... definitions
!
      dc    = ( d1 + d2 ) * 0.5D0
      ind1 = 1.0D0 / d1
      ind2 = 1.0D0 / d2
      indc = 1.0D0 / dc 

      u1n  = ( ucn + umn ) * 0.5D0

      t2n  = tg( n2 )
      t1n  = tg( n1 )
      t0n  = tg( n0 )
      tpn  = ( t1n + t0n ) * 0.5D0
      tcn  = ( t1n + t2n ) * 0.5D0
      epcn = ( ep( n1 ) + ep( n2 ) ) * 0.5D0
!
      u2n = umn - ( ucn - umn ) * ind1 * 0.5 * d2
!
! ... OUTFLOW ...

      IF( ucn < 0.D0 ) THEN
!
! ... Extrapolations
        umnn = u1n
        t1nn = tpn
! 
! ... calculation of the gas volumetric fraction at time (n+1)
! ... from the mass balance equation of solids in cell (ijk)

        epc = 0.D0
        DO is = 1, nsolid
         eps = rlk(n1,is) * inrl(is) -                                 &
             dt*inrl(is)*ind1 * (uscn(is)*rlk(n0,is) - usmn(is)*rlk(n1,is))
         epc = epc + eps
        END DO
        ep1nn = 1.D0 - epc
!
! ... Mixture Density at time (n)dt
        rm2n=0.D0
        rm1n=0.D0
        rm0n=0.D0
        DO is=1,nsolid
          rm2n=rm2n+rlk(n2,is)
          rm1n=rm1n+rlk(n1,is)
          rm0n=rm0n+rlk(n0,is)
        END DO
        rm1knn = rm1n
        rm2n = rm2n + rgp(n2)
        rm1n = rm1n + rgp(n1)
        rm0n = rm0n + rgp(n0)
!
        rmcn=(rm2n+rm1n)*0.5D0
        rmpn=(rm1n+rm0n)*0.5D0
!
        rmcnn = rmcn - dt*indc * (u1n*rmpn-u2n*rmcn)
!
! ... calculation of fluid pressure 
! ... from a mass balance equation for the mixture.
! ... Use gas velocity for mixture.
!
        mg=0.D0
        DO ig=1,ngas
          mg = mg + xgc(n1,ig) * gmw(gas_type(ig))
        END DO
        rm1nn = ep1nn*mg/(rgas*t1nn)
        p1nn = (1.D0/rm1nn) * &
               (-rm1knn+rm1n-dt*ind1*(ucn*rm0n-rm1n*umn))
        IF (epc < 1.0D-8) p1nn = p(n1)
!
! ... Calculation of the advanced-time fluid pressure from Momentum balance
! ... equation of the mixture

        p(n2) = -rmcnn*umnn + rmcn*umn - dt*indc *     &
                 ( u1n*ucn*rmpn - u2n*umn*rmcn) + dt*indc * p1nn
        p(n2)=p(n2)/(dt*indc)
!
! ... Correct non-physical pressure
        IF (p(n2) <= 0.0D0) p(n2) = p(n1)

        ep(n2) = ep(n1)
        tg(n2) = tg(n1)
        DO ig=1,ngas
          ygc(n2,ig) = ygc(n2,ig) - &
                         ucn * dt*ind2 * (ygc(n1,ig)-ygc(n2,ig))
        END DO
!
! ... extrapolation of the temperature and solid fraction to time (n+1)dt

        DO is=1,nsolid
          ! IF(nfllt.EQ.4) ws(n2,is)=ws(n1,is)
          rlk(n2,is)=rlk(n1,is)
          usmn(is)=uscn(is)
        END DO

      ELSE IF( ucn > 0.D0 ) THEN

! ... INFLOW ...
!
        zrif = zb(k) + 0.5D0 * ( dz(1) - dz(k) )  ! DOMANDA percheÃ dz(1)
        prif = p_atm(k)
        trif = t_atm(k)
        
        !rhorif = prif*gmw(6)/(rgas*trif)
        !cost = prif/(rhorif**gammaair)
        !costc = (gammaair*cost**(1.D0/gammaair))/(gammaair-1.D0)

        rhorif = rgas * trif/( prif*gmw(6))
        cost = prif**(1.D0/gammaair) * rhorif
        costc = (gammaair-1.D0)/ (gammaair*cost)
!
! ... Adiabatic inflow

        ! p(n2) = ( prif**gamn - ( u2n**2 ) / ( 2.D0 * costc ) ) ** ( 1.D0 / gamn )

        pf1   = prif**gamn
        !pf2   = ( u2n**2 ) / ( 2.D0 * costc )
        pf2   = u2n*u2n * 0.5D0 * costc 
        pfd   = ( pf1 - pf2 ) 
        IF( pfd >= 0 ) THEN
          p(n2) = pfd ** ( 1.d0 / gamn )
        ELSE
          p(n2) = prif
        END IF
!
! ... Correct non-physical pressure

        IF ( p(n2) < 0.0D0 ) p(n2) = prif
!
        ep(n2) = 1.D0
        tg(n2) = trif
        DO ig=1,ngas
          ygc(n2,ig) = ygc(n1,ig)
        END DO

      ENDIF
!
! .... Set primary variables
!
      DO is=1,nsolid
        IF( uscn(is) >= 0.D0 ) THEN
          rlk(n2,is)=rlk(n1,is)
          usmn(is)=uscn(is)
        ELSE
          rlk(n2,is)=0.0D0
          usmn(is) =0.0D0
        ENDIF
      END DO
!
      sieg(n2) = sieg(n1)
      DO is=1,nsolid
        sies(n2,is) = sies(n1,is)
      END DO
!
      RETURN
      END SUBROUTINE wsb_inout4
!-----------------------------------------------------------------------
      SUBROUTINE ent_inout6(ucn, upn, uscn, uspn, d1, d2, k)
! ... This routine computes continuous inoutflow conditions.
! ... This procedure is suited for low-Mach number regimes
! ... (incompressible flow)
! ... Works on East, North, Top boundaries

      USE atmospheric_conditions, ONLY: p_atm, t_atm, atm_ygc
      USE gas_constants, ONLY: gmw, gammaair, gamn, rgas
      USE gas_constants, ONLY: gas_type
      USE eos_gas, ONLY: mole, thermal_eosg

      REAL*8, INTENT(IN) :: ucn
      REAL*8, INTENT(INOUT) :: upn
      REAL*8, INTENT(IN) :: uscn(:)
      REAL*8, INTENT(INOUT) :: uspn(:)
      REAL*8, INTENT(IN) :: d1, d2
      INTEGER, INTENT(IN) :: k

      REAL*8 :: kfl, deltau, rls, mg
      REAL*8 :: zrif, prif, trif, rhorif
      INTEGER :: is, ig
      REAL*8 :: xgcn2(max_ngas)
!
! ... This value can be "tuned"  (0.0 <= kfl <= 1.0)
! ... 'kfl=1' corresponds to the zero-gradient condition

      kfl = 1.0D0

! ... Transport normal velocity in boundary cells
!
      deltau = upn - ucn
      upn = upn - kfl * deltau
      DO is = 1, nsolid
        deltau = uspn(is) - uscn(is)
        uspn(is) = uspn(is) - kfl * deltau
      END DO

! ... Transport particles
!
      rls = 0.D0
      DO is=1, nsolid
        IF (uscn(is) > 0.D0) THEN
          rlk(n2,is) = rlk(n2,is) - dt/d2 * uscn(is) * &
                     ( rlk(n2,is) - rlk(n1,is) )
        ELSE IF (uscn(is) < 0.D0) THEN
          rlk(n2,:) = 0.D0
        ELSE
        END IF
        rls = rls + rlk(n2,is)*inrl(is)
      END DO
      IF (rls >= 1.D0) &
          CALL error('bdry','control transport on boundary',1)

! ... Compute void fraction 
!
      ep(n2) = 1.D0 - rls

! ... Transport temperature 
!
      zrif = zb(k) + 0.5D0 * ( dz(1) - dz(k) )  ! DOMANDA perche dz(1)
      prif = p_atm(k)
      trif = t_atm(k)

      IF (ucn > 0.D0) THEN
        tg(n2) = tg(n1)
      ELSE IF (ucn < 0.D0) THEN
        tg(n2) = trif
      ELSE
      END IF
      !
      DO is=1,nsolid
        IF (uscn(is) > 0.D0) THEN
          ts(n2,is) = ts(n1,is)
        ELSE IF (uscn(is) < 0.D0) THEN
          ts(n2,is) = trif
        ELSE
        END IF
      END DO

! ... Transport gas components
!
      IF (ucn > 0.D0) THEN
        DO ig=1, ngas
          ygc(n2,ig) = ygc(n2,ig) - dt/d2 * ucn * ( ygc(n2,ig) - ygc(n1,ig) )
        END DO
      ELSE IF (ucn < 0.D0) THEN
        DO ig = 1, ngas
          ygc(n2,ig) = atm_ygc(gas_type(ig))
        END DO
      ELSE
      END IF
      CALL mole(xgcn2(:), ygc(n2,:))
      CALL thermal_eosg(rhorif, trif, prif, xgcn2(:))

! ... Transport gas bulk density
!
      IF (ucn > 0.D0) THEN
        rgp(n2) = rgp(n2) - dt/d2 * ucn * &
                     ( rgp(n2) - rgp(n1) )
      ELSE IF (ucn < 0.D0) THEN
        rgp(n2) = rhorif
      ELSE
      END IF
      rog(n2) = rgp(n2) / ep(n2)

! ... Compute new pressure by using equation of state
!
      mg = 0.D0
      DO ig = 1, ngas
        mg = mg + xgcn2(ig) * gmw(gas_type(ig))
      END DO

      p(n2) = rog(n2) * tg(n2) * ( rgas / mg )

      RETURN
      END SUBROUTINE ent_inout6
!-----------------------------------------------------------------------
      SUBROUTINE wsb_inout6(ucn, umn, uscn, usmn, d1, d2, k)
! ... This routine computes continuous inoutflow conditions.
! ... This procedure is suited for low-Mach number regimes
! ... (incompressible flow)
! ... Works on West, South, Bottom boundaries

      USE atmospheric_conditions, ONLY: p_atm, t_atm, atm_ygc
      USE gas_constants, ONLY: gmw, gammaair, gamn, rgas
      USE gas_constants, ONLY: gas_type
      USE eos_gas, ONLY: mole, thermal_eosg

      REAL*8, INTENT(IN) :: ucn
      REAL*8, INTENT(INOUT) :: umn
      REAL*8, INTENT(IN) :: uscn(:)
      REAL*8, INTENT(INOUT) :: usmn(:)
      REAL*8, INTENT(IN) :: d1, d2
      INTEGER, INTENT(IN) :: k

      REAL*8 :: kfl, deltau, rls, mg
      REAL*8 :: zrif, prif, trif, rhorif
      INTEGER :: is, ig
      REAL*8 :: xgcn2(max_ngas)
!
! ... This value can be "tuned" (0.0 <= kfl <= 1.0)
! ... (analogous to CFL number)
! ... 'kfl=1' corresponds to the zero-gradient condition
!
      kfl = 1.0D0

! ... Transport normal velocity into boundary cells
!
      deltau = ucn - umn
      umn = umn + kfl * deltau
      DO is = 1, nsolid
        deltau = uscn(is) - usmn(is)
        usmn(is) = usmn(is) - kfl * deltau
      END DO

! ... Transport particles (no g-p interaction)
!
      rls = 0.D0
      DO is=1, nsolid
        IF (usmn(is) < 0.D0) THEN
          rlk(n2,is) = rlk(n2,is) & 
                        - dt/d2 * usmn(is) * ( rlk(n1,is) - rlk(n2,is) )
        ELSE IF (usmn(is) > 0.D0) THEN
          rlk(n2,is) = 0.D0
        END IF
        rls = rls + rlk(n2,is)*inrl(is)
      END DO
      IF (rls >= 1.D0) &
          CALL error('bdry','control transport on boundary',1)

! ... Compute void fraction 
!
      ep(n2) = 1.D0 - rls

! ... Transport temperature
!
      zrif = zb(k) + 0.5D0 * ( dz(1) - dz(k) )  ! DOMANDA perche dz(1)
      prif = p_atm(k)
      trif = t_atm(k)

      IF (umn < 0.D0) THEN
        tg(n2) = tg(n1)
      ELSE IF (umn > 0.D0) THEN
        tg(n2) = trif
      ELSE
      END IF
      !
      DO is=1,nsolid
        IF (usmn(is) < 0.D0) THEN
          ts(n2,is) = ts(n1,is)
        ELSE IF (usmn(is) > 0.D0) THEN
          ts(n2,is) = trif
        ELSE
        END IF
      END DO

! ... Transport gas components
!
      IF (umn < 0.D0) THEN
        DO ig=1, ngas
          ygc(n2,ig) = ygc(n2,ig) - dt/d2 * umn * ( ygc(n1,ig) - ygc(n2,ig) )
        END DO
      ELSE IF (umn > 0.D0) THEN
        DO ig = 1, ngas
          ygc(n2,ig) = atm_ygc(gas_type(ig))
        END DO
      ELSE
      END IF
      CALL mole(xgcn2(:), ygc(n2,:))
      CALL thermal_eosg(rhorif, trif, prif, xgcn2(:))

! ... Transport gas bulk density
!
      IF (umn < 0.D0) THEN
        rgp(n2) = rgp(n2) - dt/d2 * umn * (rgp(n1) - rgp(n2))
      ELSE IF (umn > 0.D0) THEN
        rgp(n2) = rhorif
      ELSE
      END IF
      rog(n2) = rgp(n2) / ep(n2)

! ... Compute new pressure by using equation of state
!
      mg = 0.D0
      DO ig = 1, ngas
        mg = mg + xgcn2(ig) * gmw(gas_type(ig))
      END DO
!
      p(n2) = rog(n2) * tg(n2) * ( rgas / mg )
!
      RETURN
      END SUBROUTINE wsb_inout6
!-----------------------------------------------------------------------
!     I N T E R P O L A T I O N     P R O C E D U R E S 
!-----------------------------------------------------------------------
      REAL*8 FUNCTION velint(fpt, vel, ijk, cx, cy, cz)
!
! ... Interpolate velocities on a forcing point to get no-slip
! ... conditions on a solid immersed boundary 

      USE set_indexes, ONLY: ipjk, imjk, ippjk, immjk, ijpk, ipjpk,    &
        imjpk, ijmk, ipjmk, imjmk, ijppk, ijmmk, ijkp, ipjkp, imjkp,   &
        ijpkp, ijmkp, ijkm, ipjkm, imjkm, ijpkm, ijmkm, ijkpp, ijkmm
      USE immersed_boundaries, ONLY: forcing_point
      IMPLICIT NONE

      TYPE(forcing_point), INTENT(IN) :: fpt
      REAL*8, DIMENSION(:), INTENT(IN) :: vel
      REAL*8, DIMENSION(:), INTENT(IN) :: cx, cy, cz
      INTEGER, INTENT(IN) :: ijk

      INTEGER :: i, j, k
      REAL*8 :: nsx, nsy, nsz
      INTEGER :: immjkpp, ippjkpp
      INTEGER :: interp

      REAL*8 :: h          !distance between the (i,j,k)-node and boundary
      REAL*8 :: zA         !distance between the boundary and the first 
                           !external node
      REAL*8 :: zB         !distance between the boundary and the second 
                           !external node
      REAL*8 :: alpha,beta !coefficients for bilinear interpolation

      immjkpp = imjkp
      ippjkpp = ipjkp
      
      interp = fpt%int
      i      = fpt%i
      j      = fpt%j
      k      = fpt%k
      nsx    = fpt%nsl%x
      nsy    = fpt%nsl%y
      nsz    = fpt%nsl%z

      SELECT CASE (interp)

      CASE (-3)

!===========================================
!====   interpolazione lineare con i	====
!====   valori nei nodi NW-NNWW		====
!===========================================	
   
         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i-1)-nsx)**2+(cz(k+1)-nsz)**2)
         IF (h <= zA) THEN
            velint=-h/zA*vel(imjkp)
         ELSE
            zB=SQRT((cx(i-2)-nsx)**2+(cz(k+2)-nsz)**2)
            velint=-((zB-h)*vel(imjkp)+(h-zA)*vel(immjkpp))/(zB-zA)
         ENDIF

      CASE (-2)
   
!===========================================
!====   interpolazione lineare sin. con ====
!====   i valori nei nodi W-WW		====
!===========================================

         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i-1)-nsx)**2+(cz(k)-nsz)**2)
         IF (h <= zA) THEN
            velint=-h/zA*vel(imjk)
         ELSE
            zB=SQRT((cx(i-2)-nsx)**2+(cz(k)-nsz)**2)
            velint=-((zB-h)*vel(imjk)+(h-zA)*vel(immjk))/(zB-zA)
         ENDIF

      CASE (-1)
   
!===========================================
!====   interpolazione bilin. con i	====
!====   valori nei nodi W-NW-N		====
!===========================================

         alpha=(cx(i-1)-nsx)/(cx(i-1)-cx(i))
         beta=(cz(k+1)-nsz)/(cz(k+1)-cz(k))
         velint=-(alpha*(1-beta)*vel(ijkp)+(1-alpha)*(1-beta)*vel(imjkp)+ &
                 (1-alpha)*beta*vel(imjk))/(alpha*beta)
!
      CASE (0)
   
!===========================================
!====   interpolazione lineare con i	====
!====   valori nei nodi N-NN		====
!===========================================
   
         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i)-nsx)**2+(cz(k+1)-nsz)**2)
         IF (h <= zA) THEN
            velint=-h/zA*vel(ijkp)
         ELSE
            zB=SQRT((cx(i)-nsx)**2+(cz(k+2)-nsz)**2)
            velint=-((zB-h)*vel(ijkp)+(h-zA)*vel(ijkpp))/(zB-zA)

! ... Quadratic interpolation
!            velint=- ((vel(ijkp)*zB-vel(ijkpp)*zA)/(zA*zB*(zA-zB))*h**2   &
!                    +(vel(ijkpp)*zA**2-vel(ijkp)*zB**2)/(zA*zB*(zA-zB))*h)
         ENDIF

      CASE (1)
   
!===========================================
!====   interpolazione bilin. con i	====
!====   valori nei nodi E-NE-N		====
!===========================================
   
         alpha=(cx(i+1)-nsx)/(cx(i+1)-cx(i))
         beta=(cz(k+1)-nsz)/(cz(k+1)-cz(k))
         velint=-(alpha*(1-beta)*vel(ijkp)+(1-alpha)*(1-beta)*vel(ipjkp)+ &
                 (1-alpha)*beta*vel(ipjk))/(alpha*beta)
         
      CASE (2)
   
!===========================================
!====   interpolazione lineare con i	====
!====   valori nei nodi E-EE		====
!===========================================
   
         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i+1)-nsx)**2+(cz(k)-nsz)**2)
         IF (h <= zA) THEN
            velint=-h/zA*vel(ipjk)
         ELSE
            zB=SQRT((cx(i+2)-nsx)**2+(cz(k)-nsz)**2)
            velint=-((zB-h)*vel(ipjk)+(h-zA)*vel(ippjk))/(zB-zA)
         ENDIF
         
      CASE (3)
   
!===========================================
!====   interpolazione lineare con i	====
!====   valori nei nodi NE-NNEE		====
!===========================================
   
         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i+1)-nsx)**2+(cz(k+1)-nsz)**2)
         IF (h <= zA) THEN
            velint=-h/zA*vel(ipjkp)
         ELSE
            zB=SQRT((cx(i+2)-nsx)**2+(cz(k+2)-nsz)**2)
            velint=-((zB-h)*vel(ipjkp)+(h-zA)*vel(ippjkpp))/(zB-zA)
         ENDIF
   
      CASE (17)
   
!===========================================
!====   interpolazione lineare      	====
!====   esterna             		====
!===========================================
   
         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i)-nsx)**2+(cz(k+1)-nsz)**2)
         velint= + h/zA*vel(ijkp)

! ... Quadratic interpolation
!         zB=SQRT((cx(i)-nsx)**2+(cz(k+2)-nsz)**2)
!         velint= (vel(ijkp)*zB-vel(ijkpp)*zA)/(zA*zB*(zA-zB))*h**2 &
!                +(vel(ijkpp)*zA**2-vel(ijkp)*zB**2)/(zA*zB*(zA-zB))*h

      CASE DEFAULT
   
!===========================================
!====   nessuna interpolazione  ============
!===========================================
   
         velint=vel(ijk)

      END SELECT
        
      END FUNCTION velint
!----------------------------------------------------------------------
      REAL*8 FUNCTION velint3d(fpt, vel, ijk, cx, cy, cz)
!
! ... Interpolate velocities on a forcing point to get no-slip
! ... conditions on a solid immersed boundary 

      USE set_indexes
      USE immersed_boundaries, ONLY: forcing_point
      IMPLICIT NONE

      TYPE(forcing_point), INTENT(IN) :: fpt
      REAL*8, DIMENSION(:), INTENT(IN) :: vel
      REAL*8, DIMENSION(:), INTENT(IN) :: cx, cy, cz
      INTEGER, INTENT(IN) :: ijk

      INTEGER :: i, j, k
      REAL*8 :: nsx, nsy, nsz
      INTEGER :: interp, delta_i, delta_j, delta_k
      INTEGER :: index_q, index_qq
      LOGICAL :: diagonal

      REAL*8 :: h          !distance between the (i,j,k)-node and boundary
      REAL*8 :: zA         !distance between the boundary and the first 
                           !external node
      REAL*8 :: zB         !distance between the boundary and the second 
                           !external node
      REAL*8 :: alpha,beta !coefficients for bilinear interpolation

      velint3d = 0.D0

      interp = fpt%int
      i      = fpt%i
      j      = fpt%j
      k      = fpt%k
      nsx    = fpt%nsl%x
      nsy    = fpt%nsl%y
      nsz    = fpt%nsl%z
!
! ... LINEAR interpolation criteria: 
! ... 0: top
! ... 1: west
! ... 2: south-west
! ... 3: south
! ... 4: south-east
! ... 5: east
! ... 6: north-east
! ... 7: north
! ... 8: north-west
! ... 17: external top

      SELECT CASE (interp)

        CASE (0)
        delta_i = 0
        delta_j = 0
        delta_k = 1
        index_q  = ijkp
        index_qq = ijkpp

        CASE (1)
        delta_i = -1
        delta_j = 0
        delta_k = 0
        index_q  = imjk
        index_qq = immjk

        CASE (2)
        delta_i = -1
        delta_j = -1
        delta_k = 0
        index_q  = imjmk

        CASE (3)
        delta_i = 0
        delta_j = -1
        delta_k = 0
        index_q  = ijmk
        index_qq = ijmmk

        CASE (4)
        delta_i = 1
        delta_j = -1
        delta_k = 0
        index_q  = ipjmk

        CASE (5)
        delta_i = 1
        delta_j = 0
        delta_k = 0
        index_q  = ipjk
        index_qq = ippjk

        CASE (6)
        delta_i = 1
        delta_j = 1
        delta_k = 0
        index_q  = ipjpk

        CASE (7)
        delta_i = 0
        delta_j = 1
        delta_k = 0
        index_q  = ijpk
        index_qq = ijppk

        CASE (8)
        delta_i = -1
        delta_j = 1
        delta_k = 0
        index_q  = imjpk

        CASE (17)
        delta_i = 0
        delta_j = 0
        delta_k = 1
        index_q  = ijkp
        index_qq = ijkpp

        CASE DEFAULT
        delta_i = 0
        delta_j = 0
        delta_k = 0
        index_q  = ijk
        index_qq = ijk


      END SELECT

      diagonal = ( (MOD(interp,2) == 0) .AND. (interp/=0) )

      h  = SQRT( (cx(i)-nsx)**2 + (cy(j)-nsy)**2 + (cz(k)-nsz)**2 )
      
      zA = SQRT( (cx(i+delta_i)-nsx)**2 + (cy(j+delta_j)-nsy)**2 + &
                 (cz(k+delta_k)-nsz)**2 )

      IF (interp == 17) THEN
         velint3d = h/zA*vel(index_q)         
      ELSEIF (h <= zA .OR. diagonal) THEN
         velint3d = -h/zA*vel(index_q)
      ELSE 
         zB=SQRT( (cx(i+2*delta_i)-nsx)**2 + (cy(j+2*delta_j)-nsy)**2 +  &
                  (cz(k+2*delta_k)-nsz)**2 )
         velint3d = -((zB-h)*vel(index_q)+(h-zA)*vel(index_qq))/(zB-zA)
      ENDIF
        
      RETURN
      END FUNCTION velint3d
!----------------------------------------------------------------------
      END MODULE boundary_conditions
!-----------------------------------------------------------------------
