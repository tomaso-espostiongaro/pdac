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
      USE interpolate_fields, ONLY: velint, velint3d
      USE eos_gas, ONLY: ygc
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p, ep
      USE time_parameters, ONLY: dt, time, sweep
!
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
      USE io_files, ONLY: tempunit, testunit
      USE domain_mapping, ONLY: ncint, myijk, meshinds
      USE grid, ONLY: flag, x, y, z, xb, yb, zb
      USE grid, ONLY: slip_wall, noslip_wall
      USE grid, ONLY: immb_cell, filled_cell_1, filled_cell_2
      USE grid, ONLY: free_io, nrfree_io
      USE grid, ONLY: inlet_cell, vent_cell, dome_cell, fluid, bl_cell
      USE immersed_boundaries, ONLY: fptx, fpty, fptz, forcing_point
      USE immersed_boundaries, ONLY: numx, numy, numz, immb
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE inflow_outflow, ONLY: n0, n1, n2
      USE inflow_outflow, ONLY: ent_inout4, wsb_inout4, extrapolate
      USE parallel, ONLY: mpime, root
      USE set_indexes, ONLY: subscr
      USE set_indexes, ONLY: ipjk, imjk, ippjk, immjk, ijpk, ipjpk,    &
        imjpk, ijmk, ipjmk, imjmk, ijppk, ijmmk, ijkp, ipjkp, imjkp,   &
        ijpkp, ijmkp, ijkm, ipjkm, imjkm, ijpkm, ijmkm, ijkpp, ijkmm
      USE time_parameters, ONLY: tau1, tau2
      USE vent_conditions, ONLY: update_vent_cell, random_switch, irand, &
                                 update_inlet_cell, grow_vent_cell, &
                                 grow_inlet_cell
!
      IMPLICIT NONE
!
      INTEGER :: ijk, i, j, k, imesh, ig, is, nph
      INTEGER :: fp, np
      REAL*8 :: d0, d1, d2 
      REAL*8 :: vel(max_nsolid+1)
      INTEGER :: fx, fy, fz
      INTEGER :: nfptx, nfpty, nfptz
      LOGICAL :: forced
      REAL*8, ALLOCATABLE :: tfptx(:), tfpty(:), tfptz(:)
      INTEGER :: indexq
!
      vel = 0.D0
      !
      ! ... Compute a random factor for vent antialiasing
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
            IF (flag(ijk) == vent_cell) THEN
                    ! ... Vent Antialiasing
                    CALL update_vent_cell(ijk,imesh,sweep)
            END IF
          END IF
!
          IF (tau1 > 0.D0 .OR. tau2 > 0.D0) THEN
            IF (flag(ijk) == inlet_cell) THEN
                  CALL grow_inlet_cell(ijk,imesh,sweep)
            ELSE IF (flag(ijk) == vent_cell) THEN
                  CALL grow_vent_cell(ijk,imesh,sweep)
            END IF
          END IF
!
! ... In the immersed boundaries, update the velocity through linear
! ... interpolations 
!
        IF( flag(ijk) == immb_cell .OR. flag(ijk) == filled_cell_1 .OR. &
              flag(ijk) == filled_cell_2) THEN
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
              vel(:) = velint(fptx(fx), ug, us, ijk, xb, y, z)  
              DO nph = 1, nsolid+1
                fptx(fx)%vel(nph) = vel(nph)
              END DO
              !
              ! ... Set the homogeneous Neumann conditions in non-resolved forcing points
              ! ... (zero-gradient)
              IF( fptx(fx)%int >= 20 ) CALL hneumann(ijk,ipjk)
              !
              ! ... Initialize x-velocity in the forced points
              ug(ijk) = vel(1)
              DO is=1,nsolid
                us(ijk,is) = vel(1+is)
              END DO
            END IF

            IF( fz/=0 ) THEN
              vel(:) = velint(fptz(fz), wg, ws, ijk, x, y, zb)  
              DO nph = 1, nsolid+1
                fptz(fz)%vel(nph) = vel(nph)
              END DO
              !
              ! ... Set the homogeneous Neumann conditions in non-resolved forcing points
              IF( fptz(fz)%int >= 20 ) CALL hneumann(ijk,ijkp)
              !
              ! ... Initialize z-velocity in the forced points
              wg(ijk) = vel(1)
              DO is=1,nsolid
                ws(ijk,is) = vel(1+is)
              END DO
            END IF

          ELSE IF (job_type == '3D') THEN

            IF( fx/=0 ) THEN
              vel(:) = velint3d(fptx(fx), ug, us, ijk, xb, y, z, indexq)  
              DO nph = 1, nsolid+1
                fptx(fx)%vel(nph) = vel(nph)
              END DO
              !
              ! ... Set the homogeneous Neumann conditions in non-resolved forcing points
              IF( fptx(fx)%int >= 20 ) CALL hneumann(ijk,ipjk)
              !
              ! ... Initialize x-velocity in the forced points
              ug(ijk) = vel(1)
              DO is=1,nsolid
                us(ijk,is) = vel(1+is) 
              END DO
            END IF
            
            IF( fy/=0 ) THEN
              vel(:) = velint3d(fpty(fy), vg, vs, ijk, x, yb, z, indexq)  
              DO nph = 1, nsolid+1
                fpty(fy)%vel(nph) = vel(nph)
              END DO

              ! ... Set the homogeneous Neumann conditions in non-resolved forcing points
              IF( fpty(fy)%int >= 20 ) CALL hneumann(ijk,ijpk)
              !
              ! ... Initialize y-velocity in the forced points
              vg(ijk) = vel(1)
              DO is=1,nsolid
                vs(ijk,is) = vel(1+is)
              END DO
            END IF
            
            IF( fz/=0 ) THEN
              vel(:) = velint3d(fptz(fz), wg, ws, ijk, x, y, zb, indexq)  
              DO nph = 1, nsolid+1
                fptz(fz)%vel(nph) = vel(nph)
              END DO

              ! ... Set the homogeneous Neumann conditions in non-resolved forcing points
              IF( fptz(fz)%int >= 20 ) CALL hneumann(ijk,ijkp)
              !
              ! ... Initialize z-velocity in the forced points
              wg(ijk) = vel(1)
              DO is=1,nsolid
                ws(ijk,is) = vel(1+is)
              END DO
            END IF
          END IF
!
! ... In fluid cells update the neighbours on boundaries
!
       ELSE IF( flag(ijk)==fluid .OR. flag(ijk)==dome_cell .OR. flag(ijk)==bl_cell) THEN
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
              d0 = dx(i-1)
              d1 = dx(i)
	      d2 = dx(i+1)

              ! ... Compute the normal component of the velocities 
              ! ... and scalar fields
              !
              CALL extrapolate( ug(n0), ug(n1), ug(n2), us(n0,:), us(n1,:), us(n2,:), d0, d1, d2, k)

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
              d0 = dx(i+1)
              d1 = dx(i)
              d2 = dx(i-1)

              ! ... Compute the normal component of the velocities 
              ! ... and scalar fields
              !
              CALL extrapolate( ug(n0), ug(n1), ug(n2), us(n0,:), us(n1,:), us(n2,:), d0, d1, d2, k)

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
                CALL extrapolate( vg(n0), vg(n1), vg(n2), vs(n0,:), vs(n1,:), vs(n2,:), d0, d1, d2, k)

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
                CALL extrapolate( vg(n0), vg(n1), vg(n2), vs(n0,:), vs(n1,:), vs(n2,:), d0, d1, d2, k)

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
              CALL extrapolate( wg(n0), wg(n1), wg(n2), ws(n0,:), ws(n1,:), ws(n2,:), d0, d1, d2, k)
                             
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
      IF (lpr > 2 .AND. immb >= 1) THEN
        !
        nfptx = SIZE(fptx)
        ALLOCATE(tfptx(nfptx))
        tfptx(:) = fptx(:)%vel(1)
        CALL parallel_sum_real(tfptx(:),nfptx)
        !
        IF (job_type == '3D') THEN
          nfpty = SIZE(fpty)
          ALLOCATE(tfpty(nfpty))
          tfpty(:) = fpty(:)%vel(1)
          CALL parallel_sum_real(tfpty(:),nfpty)
        END IF
        !
        nfptz = SIZE(fptz)
        ALLOCATE(tfptz(nfptz))
        tfptz(:) = fptz(:)%vel(1)
        CALL parallel_sum_real(tfptz(:),nfptz)
        !
        IF (mpime == root) THEN
          OPEN(UNIT=tempunit,FILE='fptx.dat',STATUS='UNKNOWN')
          DO np = 1, SIZE(fptx)
            WRITE(tempunit,32) np, fptx(np)%i, fptx(np)%j, fptx(np)%k, fptx(np)%int, fptx(np)%nsl, tfptx(np)
          END DO
          CLOSE(tempunit)
          IF (job_type == '3D') THEN
            OPEN(UNIT=tempunit,FILE='fpty.dat',STATUS='UNKNOWN')
            DO np = 1, SIZE(fpty)
              WRITE(tempunit,32) np, fpty(np)%i, fpty(np)%j, fpty(np)%k, fpty(np)%int, fpty(np)%nsl, tfpty(np)
            END DO
            CLOSE(tempunit)
          END IF
          OPEN(UNIT=tempunit,FILE='fptz.dat',STATUS='UNKNOWN')
          DO np = 1, SIZE(fptz)
            WRITE(tempunit,32) np, fptz(np)%i, fptz(np)%j, fptz(np)%k, fptz(np)%int, fptz(np)%nsl, tfptz(np)
          END DO
          CLOSE(tempunit)
          !
        END IF
  32    FORMAT(5(I6),4(F12.3))
        !
        DEALLOCATE(tfptx)
        IF (job_type == '3D') DEALLOCATE(tfpty)
        DEALLOCATE(tfptz)
      END IF
!
      RETURN
      END SUBROUTINE boundary
!----------------------------------------------------------------------
      SUBROUTINE hneumann(m,n)
      USE grid, ONLY: flag, filled_cell_1, filled_cell_2, fc2
      USE pressure_epsilon, ONLY: p, ep
      USE eos_gas, ONLY: ygc
      USE gas_solid_density, ONLY: rlk, rgp
      USE gas_solid_temperature, ONLY: tg, ts, sieg, sies
      USE io_files, ONLY: testunit
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: m,n
      LOGICAL :: fc2_
!
! ... m is the local cell index; 
! ... n is the index of the first neighbour used for the external forcing
!
      IF (flag(m) == filled_cell_1) THEN
        p(m) = p(n)
        ep(m) = ep(n)
        rlk(m,:) = rlk(n,:)
        rgp(m) = rgp(n)
        tg(m) = tg(n)
        ts(m,:) = ts(n,:)
        sieg(m) = sieg(n)
        sies(m,:) = sies(n,:)
        ygc(m,:) = ygc(n,:)
      ELSE IF (flag(n) == filled_cell_1) THEN
        p(n) = p(m)
        ep(n) = ep(m)
        rlk(n,:) = rlk(m,:)
        rgp(n) = rgp(m)
        tg(n) = tg(m)
        ts(n,:) = ts(m,:)
        sieg(n) = sieg(m)
        sies(n,:) = sies(m,:)
        ygc(n,:) = ygc(m,:)
      END IF
!      
      IF (flag(m) == filled_cell_2) THEN
        fc2_ = fc2(m)
        IF (.NOT.fc2_) THEN
          p(m) = 0.5D0 * p(n)
          ep(m) = 0.5D0 * ep(n)
          rlk(m,:) = 0.5D0 * rlk(n,:)
          rgp(m) = 0.5D0 * rgp(n)
          tg(m) = 0.5D0 * tg(n)
          ts(m,:) = 0.5D0 * ts(n,:)
          sieg(m) = 0.5D0 * sieg(n)
          sies(m,:) = 0.5D0 * sies(n,:)
          ygc(m,:) = 0.5D0 * ygc(n,:)
        ELSE
          p(m) = p(m) + 0.5D0 * p(n)
          ep(m) = ep(m) + 0.5D0 * ep(n)
          rlk(m,:) = rlk(m,:) + 0.5D0 * rlk(n,:)
          rgp(m) = rgp(m) + 0.5D0 * rgp(n)
          tg(m) = tg(m) + 0.5D0 * tg(n)
          ts(m,:) = ts(m,:) + 0.5D0 * ts(n,:)
          sieg(m) = sieg(m) + 0.5D0 * sieg(n)
          sies(m,:) = sies(m,:) + 0.5D0 * sies(n,:)
          ygc(m,:) = ygc(m,:) + 0.5D0 * ygc(n,:)
        END IF
        fc2(m) = .NOT.fc2_
      ELSE IF (flag(n) == filled_cell_2) THEN
        fc2_ = fc2(n)
        IF (.NOT.fc2_) THEN
          p(n) = 0.5D0 * p(m)
          ep(n) = 0.5D0 * ep(m)
          rlk(n,:) = 0.5D0 * rlk(m,:)
          rgp(n) = 0.5D0 * rgp(m)
          tg(n) = 0.5D0 * tg(m)
          ts(n,:) = 0.5D0 * ts(m,:)
          sieg(n) = 0.5D0 * sieg(m)
          sies(n,:) = 0.5D0 * sies(m,:)
          ygc(n,:) = 0.5D0 * ygc(m,:)
        ELSE
          p(n) = p(n) + 0.5D0 * p(m)
          ep(n) = ep(n) + 0.5D0 * ep(m)
          rlk(n,:) = rlk(n,:) + 0.5D0 * rlk(m,:)
          rgp(n) = rgp(n) + 0.5D0 * rgp(m)
          tg(n) = tg(n) + 0.5D0 * tg(m)
          ts(n,:) = ts(n,:) + 0.5D0 * ts(m,:)
          sieg(n) = sieg(n) + 0.5D0 * sieg(m)
          sies(n,:) = sies(n,:) + 0.5D0 * sies(m,:)
          ygc(n,:) = ygc(n,:) + 0.5D0 * ygc(m,:)
        END IF
        fc2(n) = .NOT.fc2_
      END IF 
!
      RETURN
      END SUBROUTINE hneumann
!----------------------------------------------------------------------
      END MODULE boundary_conditions
!-----------------------------------------------------------------------
