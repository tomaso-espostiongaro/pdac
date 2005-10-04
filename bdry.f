!-----------------------------------------------------------------------
      MODULE boundary_conditions
!-----------------------------------------------------------------------
!
      USE atmospheric_conditions, ONLY: gravx, gravy, gravz
      USE dimensions
      USE gas_solid_density, ONLY: rgp, rlk, rog
      USE gas_solid_temperature, ONLY: sieg, tg, sies, ts
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE grid, ONLY: dx, dy, dz, z
      USE eos_gas, ONLY: xgc, ygc
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p, ep
      USE time_parameters, ONLY: dt, time, sweep

      IMPLICIT NONE

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
      USE domain_mapping, ONLY: ncint, myijk, meshinds
      USE grid, ONLY: flag, x, y, z, xb, yb, zb
      USE grid, ONLY: slip_wall, noslip_wall, fluid, immb_cell, filled_cell
      USE grid, ONLY: free_io, nrfree_io, inlet_cell, vent_cell
      USE immersed_boundaries, ONLY: fptx, fpty, fptz
      USE immersed_boundaries, ONLY: numx, numy, numz, immb
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE inflow_outflow, ONLY: n0, n1, n2
      USE inflow_outflow, ONLY: ent_inout4, ent_inout6, wsb_inout4, wsb_inout6
      USE interpolate_fields, ONLY: velint, velint3d
      USE parallel, ONLY: mpime, root
      USE set_indexes, ONLY: subscr
      USE set_indexes, ONLY: ipjk, imjk, ippjk, immjk, ijpk, ipjpk,    &
        imjpk, ijmk, ipjmk, imjmk, ijppk, ijmmk, ijkp, ipjkp, imjkp,   &
        ijpkp, ijmkp, ijkm, ipjkm, imjkm, ijpkm, ijmkm, ijkpp, ijkmm
      USE vent_conditions, ONLY: update_vent_cell, random_switch, irand, &
                                 update_inlet_cell
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
! ... Loop over the mesh and check boundaries
!
      mesh_loop: DO ijk = 1, ncint
        CALL subscr(ijk)
        CALL meshinds(ijk,imesh,i,j,k)

        !
        fx = 0
        fy = 0 
        fz = 0
        forced = .FALSE.
!
! ... Update inlet cells for non-stationnary boundary conditions
!
        IF (irand >= 1) THEN
          IF (flag(ijk) == inlet_cell) THEN
                  ! ... Impose a random perturbation at inlet
            CALL update_inlet_cell(ijk)
          ELSE IF (flag(ijk) == vent_cell) THEN
                  ! ... Modify randomly the vent cells laying on the 
                  ! ... vent rim to get the prescribed (averaged) mass flow-rate
            CALL update_vent_cell(ijk,imesh,sweep)
          END IF
        END IF
!
! ... In the immersed boundaries, update the velocity through linear
! ... interpolations 
!
        IF( flag(ijk) == immb_cell .OR. flag(ijk) == filled_cell) THEN
          !
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
              IF( fptx(fx)%int >= 20 ) THEN
                      IF (flag(ijk) == filled_cell) THEN
                              p(ijk) = p(ipjk)
                      ELSE
                              p(ipjk) = p(ijk)
                      END IF
              END IF
              !
              ! ... Initialize x-velocity in the forced points
              ug(ijk) = vel
              DO is=1,nsolid
                us(ijk,is) = vel
              END DO
            END IF

            IF( fz/=0 ) THEN
              vel = velint(fptz(fz), wg, ijk, x, y, zb)  
              fptz(fz)%vel = vel
              !
              ! ... Set the pressure in non-resolved forcing points
              IF( fptz(fz)%int >= 20 ) THEN
                      p(ijk) = p(ijkp)
              END IF
              !
              ! ... Initialize z-velocity in the forced points
              wg(ijk) = vel
              DO is=1,nsolid
                ws(ijk,is) = vel
              END DO
            END IF

          ELSE IF (job_type == '3D') THEN

            IF( fx/=0 ) THEN
              vel = velint3d(fptx(fx), ug, ijk, xb, y, z)  
              fptx(fx)%vel = vel
              !
              ! ... Set the pressure in non-resolved forcing points
              IF( fptx(fx)%int >= 20 ) THEN
                      IF (flag(ijk) == filled_cell) THEN
                              p(ijk) = p(ipjk)
                      ELSE
                              p(ipjk) = p(ijk)
                      END IF
              END IF
              !
              ! ... Initialize x-velocity in the forced points
              ug(ijk) = vel
              DO is=1,nsolid
                us(ijk,is) = vel
              END DO
            END IF
            
            IF( fy/=0 ) THEN
              vel = velint3d(fpty(fy), vg, ijk, x, yb, z)  
              fpty(fy)%vel = vel

              ! ... Set the pressure in non-resolved forcing points
              IF( fpty(fy)%int >= 20 ) THEN
                      IF (flag(ijk) == filled_cell) THEN
                              p(ijk) = p(ijpk)
                      ELSE
                              p(ijpk) = p(ijk)
                      END IF
              END IF
              
              ! ... Initialize y-velocity in the forced points
              vg(ijk) = vel
              DO is=1,nsolid
                vs(ijk,is) = vel
              END DO
            END IF
            
            IF( fz/=0 ) THEN
              vel = velint3d(fptz(fz), wg, ijk, x, y, zb)  
              fptz(fz)%vel = vel

              ! ... Set the pressure in non-resolved forcing points
              IF( fptz(fz)%int >= 20 ) THEN
                      p(ijk) = p(ijkp)
              END IF
              
              ! ... Initialize z-velocity in the forced points
              wg(ijk) = vel
              DO is=1,nsolid
                ws(ijk,is) = vel
              END DO
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
	      IF (job_type == '3D') THEN
                vg(n2)   = vg(n1)
                vs(n2,:) = vs(n1,:)
              END IF
              wg(n2)   = wg(n1)
              ws(n2,:) = ws(n1,:)
!
            CASE (noslip_wall)
!
              ug(n2)   = 0.D0
              us(n2,:) = 0.D0
              
              IF(.NOT.BTEST(flag(ipjkp),0)) THEN
              !IF(flag(ipjkp)/=fluid) THEN
                wg(n2)   = -wg(n1)
                DO is = 1, nsolid
                  IF (rlk(ijk,is) > 0.D0) ws(n2,is) = -ws(n1,is)
                END DO
              END IF

              IF (.NOT.BTEST(flag(ipjpk),0) .AND. job_type == '3D') THEN
              !IF (flag(ipjpk) /= fluid .AND. job_type == '3D') THEN
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
	      IF (job_type == '3D') THEN
                vg(n2)   = vg(n1)
                vs(n2,:) = vs(n1,:)
	      END IF
              wg(n2)   = wg(n1)
              ws(n2,:) = ws(n1,:)
!
            CASE (noslip_wall)

              ug(n2)   = 0.D0
              us(n2,:) = 0.D0

              IF (.NOT.BTEST(flag(imjkp),0) ) THEN
              !IF ( flag(imjkp) /= fluid ) THEN
                wg(n2)   = -wg(n1)
                DO is = 1, nsolid
                  IF (rlk(ijk,is) > 0.D0) ws(n2,is) = -ws(n1,is)
                END DO
	      END IF

              IF ( .NOT.BTEST(flag(imjkp),0) .AND. job_type == '3D') THEN
              !IF ( flag(imjkp) /= fluid .AND. job_type == '3D') THEN
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
                
                IF(.NOT.BTEST(flag(ipjpk),0)) THEN
                !IF(flag(ipjpk) /= fluid) THEN
                  ug(n2)   = -ug(n1)
                  DO is = 1, nsolid
                    IF (rlk(ijk,is) > 0.D0) us(n2,is) = -us(n1,is)
                  END DO
                END IF

                IF(.NOT.BTEST(flag(ijpkp),0)) THEN
                !IF(flag(ijpkp) /= fluid) THEN
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

                IF(.NOT.BTEST(flag(ipjmk),0)) THEN
                !IF(flag(ipjmk) /= fluid) THEN
                  ug( n2 ) = -ug( n1 )
                  DO is = 1, nsolid
                    IF (rlk(ijk,is) > 0.D0) us(n2,is) = -us(n1,is)
                  END DO
                  
                END IF
                IF(.NOT.BTEST(flag(ijmkp),0)) THEN
                !IF(flag(ijmkp) /= fluid) THEN
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

              IF (.NOT.BTEST(flag(ipjkp),0)) THEN
              !IF (flag(ipjkp) /= fluid) THEN
                ug(n2)   = -ug(n1)
                DO is = 1, nsolid
                  IF (rlk(ijk,is) > 0.D0) us(n2,is) = -us(n1,is)
                END DO
              END IF

	      IF (.NOT.BTEST(flag(ijpkp),0) .AND. job_type == '3D') THEN
	      !IF (flag(ijpkp) /= fluid .AND. job_type == '3D') THEN
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
      END MODULE boundary_conditions
!-----------------------------------------------------------------------
