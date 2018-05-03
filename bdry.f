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
      USE eos_gas, ONLY: ygc
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p, ep
      USE time_parameters, ONLY: dt, time, sweep
!
      IMPLICIT NONE
      INTEGER :: ord
      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE boundary
!
! ... This routine computes (x,y,z) boundary conditions 
!
      USE control_flags, ONLY: job_type, lpr
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE io_files, ONLY: tempunit, testunit, logunit
      USE domain_mapping, ONLY: ncint, myijk, meshinds
      USE grid, ONLY: flag, x, y, z, xb, yb, zb
      USE grid, ONLY: slip_wall, noslip_wall
      USE grid, ONLY: immb_cell, filled_cell_1, filled_cell_2
      USE grid, ONLY: free_io, zero_grad
      USE grid, ONLY: inlet_cell, vent_cell, dome_cell, fluid, bl_cell
      USE immersed_boundaries, ONLY: fptx, fpty, fptz, forcing_point
      USE immersed_boundaries, ONLY: numx, numy, numz, immb, faces
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE inflow_outflow, ONLY: n0, n1, n2
      USE inflow_outflow, ONLY: inoutflow, zerograd
      USE interpolate_fields, ONLY: velint, velint3d, extrapolate, hn
      USE mass_sink, ONLY: sink_update, sink, isink
      USE parallel, ONLY: mpime, root
      USE set_indexes, ONLY: subscr
      USE set_indexes, ONLY: ipjk, imjk, ippjk, immjk, ijpk, ipjpk,    &
        imjpk, ijmk, ipjmk, imjmk, ijppk, ijmmk, ijkp, ipjkp, imjkp,   &
        ijpkp, ijmkp, ijkm, ipjkm, imjkm, ijpkm, ijmkm, ijkpp, ijkmm
      USE time_parameters, ONLY: tau1, tau2
      USE vent_conditions, ONLY: inlet_velocity_fluctuations, random_switch, irand, iali
      USE vent_conditions, ONLY: grow_inlet_cell, update_vent_cell
      USE vent_conditions, ONLY: inlet_profile, vent_index
!
      IMPLICIT NONE
!
      INTEGER :: ijk, i, j, k, imesh, ig, is, nph, n, ventn
      INTEGER :: fp, np
      REAL*8 :: d0, d1, d2 
      REAL*8 :: vel(max_nsolid+1)
      INTEGER :: fx, fy, fz
      INTEGER :: nfptx, nfpty, nfptz
      LOGICAL :: forced
      REAL*8, ALLOCATABLE :: tfptx(:), tfpty(:), tfptz(:)
      REAL*8 :: gnorm
      INTEGER :: indexq
!
      vel = 0.D0
      !
      ! ... Compute a random factor for vent antialiasing
      IF (irand >= 1) CALL random_switch(sweep)
!
! ... Initialize sink array
!
      IF (isink > 0) sink = 0.D0
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
              CALL vent_index(imesh,n)
              ! ... random antialias
              IF (iali == 5) CALL update_vent_cell(ijk,n)
              ! ... velocity fluctuations
              IF (inlet_profile >= 0) CALL inlet_velocity_fluctuations(ijk,n)
            END IF
          END IF
!
          IF (tau1 > 0.D0 .OR. tau2 > 0.D0) THEN
            IF (flag(ijk) == inlet_cell .OR. flag(ijk) == vent_cell) THEN
                  CALL grow_inlet_cell(ijk,imesh,sweep)
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
          IF (job_type == JOB_TYPE_3D) fy = numy(ijk)
          fz = numz(ijk)
          forced = (fx/=0 .OR. fy/=0 .OR. fz/=0)
          IF (.NOT.forced) &
            CALL error('boundary','control forcing points',1)
          ! 
          ! ... Compute the pseudo-velocities
          ! ... that are used in the "immersed boundary" technique ...
          !
          IF (job_type == JOB_TYPE_2D) THEN

            IF( fx/=0 ) THEN
              vel(:) = velint(fptx(fx), ug, us, ijk, xb, y, z)  
              DO nph = 1, nsolid+1
                fptx(fx)%vel(nph) = vel(nph)
              END DO
              !
              ! ... Set the homogeneous Neumann conditions in non-resolved forcing points
              ! ... (zero-gradient)
              IF( fptx(fx)%int >= 20 ) CALL hn(ijk,ipjk)
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
              IF( fptz(fz)%int >= 20 ) CALL hn(ijk,ijkp)
              !
              ! ... Initialize z-velocity in the forced points
              wg(ijk) = vel(1)
              DO is=1,nsolid
                ws(ijk,is) = vel(1+is)
              END DO
            END IF

          ELSE IF (job_type == JOB_TYPE_3D) THEN

            IF( fx/=0 ) THEN
              vel(:) = velint3d(fptx(fx), ug, us, ijk, xb, y, z, indexq)  
              DO nph = 1, nsolid+1
                fptx(fx)%vel(nph) = vel(nph)
              END DO
              !
              ! ... Set the homogeneous Neumann conditions in non-resolved forcing points
              IF( fptx(fx)%int >= 20 ) CALL hn(ijk,ipjk)
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
              IF( fpty(fy)%int >= 20 ) CALL hn(ijk,ijpk)
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
              IF( fptz(fz)%int >= 20 ) CALL hn(ijk,ijkp)
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
              ug(n1)   = 0.D0
              us(n1,:) = 0.D0
              ug(n2)   = 0.D0
              us(n2,:) = 0.D0
!
            IF (job_type == JOB_TYPE_3D) THEN
                vg(n2)   = vg(n1)
                vs(n2,:) = vs(n1,:)
              END IF
              wg(n2)   = wg(n1)
              ws(n2,:) = ws(n1,:)
!
            CASE (noslip_wall)
!
              ug(n1)   = 0.D0
              us(n1,:) = 0.D0
              ug(n2)   = 0.D0
              us(n2,:) = 0.D0
              
              IF(.NOT.BTEST(flag(ipjkp),0)) THEN
              !IF(flag(ipjkp)/=fluid) THEN
                wg(n2)   = -wg(n1)
                ws(n2,:) = -ws(n1,:)
              END IF

              IF (.NOT.BTEST(flag(ipjpk),0) .AND. job_type == JOB_TYPE_3D) THEN
              !IF (flag(ipjpk) /= fluid .AND. job_type == JOB_TYPE_3D) THEN
                vg(n2)   = -vg(n1)
                vs(n2,:) = -vs(n1,:)
              END IF
!
            CASE (free_io, zero_grad) ! Zero gradient of all quantities
              d0 = dx(i-1)
              d1 = dx(i)
              d2 = dx(i+1)

              ! ... Zero gradient normal velocities
              ug(n1) = ug(n0)
              us(n1,:) = us(n0,:)
              ug(n2)   = ug(n1)
              us(n2,:) = us(n1,:)
!
              ! ... Tangential velocities
              IF (job_type == JOB_TYPE_3D) THEN
                vg(n2)   = vg(n1)
                vs(n2,:) = vs(n1,:)
              END IF
              wg(n2)   = wg(n1)
              ws(n2,:) = ws(n1,:)
              !
              ! ... outlet rlk, rgp, rog, sieg, sies, tg, ts, ep
              gnorm = gravx * 0.5D0 * (dx(i+1) + dx(i))
              IF (flag(n2) == zero_grad) gnorm = 0.D0
              CALL inoutflow(ug(n1),us(n1,:),n2,n1,k,gnorm)

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
              ug(n2)   = 0.D0
              us(n2,:) = 0.D0
!
              IF (job_type == JOB_TYPE_3D) THEN
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
                ws(n2,:) = -ws(n1,:)
              END IF

              IF ( .NOT.BTEST(flag(imjkp),0) .AND. job_type == JOB_TYPE_3D) THEN
              !IF ( flag(imjkp) /= fluid .AND. job_type == JOB_TYPE_3D) THEN
                vg(n2)   = -vg(n1)
                vs(n2,:) = -vs(n1,:)
            END IF
!
            CASE (free_io, zero_grad)
!	    
              d0 = dx(i+1)
              d1 = dx(i)
              d2 = dx(i-1)
              !
              ! ... Normal velocities
              ug(n2) = ug(n1)
              us(n2,:) = us(n1,:)
              !
              ! ... Tangential velocities
              IF (job_type == JOB_TYPE_3D) THEN
                vg(n2)   = vg(n1)
                vs(n2,:) = vs(n1,:)
              END IF
              wg(n2)   = wg(n1)
              ws(n2,:) = ws(n1,:)
              !
              gnorm = gravx * 0.5D0 * (dx(i-1) + dx(i))
              IF (flag(n2) == zero_grad) gnorm = 0.D0
              CALL inoutflow(ug(n2),us(n2,:),n2,n1,k,gnorm)
!              
            CASE DEFAULT
!
              CONTINUE

            END SELECT

            IF (job_type == JOB_TYPE_3D) THEN
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
                vg(n1)   = 0.D0
                vs(n1,:) = 0.D0
                vg(n2)   = 0.D0
                vs(n2,:) = 0.D0
!  
                ug(n2)   = ug(n1)
                us(n2,:) = us(n1,:)
                wg(n2)   = wg(n1)
                ws(n2,:) = ws(n1,:)
!  
              CASE (noslip_wall)
!  
                vg(n1)   = 0.D0
                vs(n1,:) = 0.D0
                vg(n2)   = 0.D0
                vs(n2,:) = 0.D0
                
                IF(.NOT.BTEST(flag(ipjpk),0)) THEN
                !IF(flag(ipjpk) /= fluid) THEN
                  ug(n2)   = -ug(n1)
                  us(n2,:) = -us(n1,:)
                END IF

                IF(.NOT.BTEST(flag(ijpkp),0)) THEN
                !IF(flag(ijpkp) /= fluid) THEN
                  wg(n2)   = -wg(n1)
                  ws(n2,:) = -ws(n1,:)
                END IF
!  
              CASE (free_io, zero_grad)
!  
                d0 = dy(j-1)
                d1 = dy(j)
                d2 = dy(j+1)

                ! ... Normal velocities
                vg(n1)   = vg(n0)
                vs(n1,:) = vs(n0,:)
                vg(n2)   = vg(n1)
                vs(n2,:) = vs(n1,:)
!
                ! ... Tangential velocities
                ug(n2)   = ug(n1)
                us(n2,:) = us(n1,:)
                wg(n2)   = wg(n1)
                ws(n2,:) = ws(n1,:)

                ! ... Outlet scalars
                gnorm = gravy * 0.5D0 * (dy(j) + dy(j+1))
                IF (flag(n2) == zero_grad) gnorm = 0.D0
                CALL inoutflow(vg(n1),vs(n1,:),n2,n1,k,gnorm)
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
                ! ... Normal velocities
                !
                vg(n2)   = 0.D0
                vs(n2,:) = 0.D0
!
                ! ... Tangential velocities
                !
                ug(n2)   = ug(n1)
                us(n2,:) = us(n1,:)
                wg(n2)   = wg(n1)
                ws(n2,:) = ws(n1,:)
!  
              CASE (noslip_wall)
!
                ! ... Normal velocities
                !
                vg(n2)   = 0.D0
                vs(n2,:) = 0.D0

                ! ... Tangential velocities
                !
                IF(.NOT.BTEST(flag(ipjmk),0)) THEN
                !IF(flag(ipjmk) /= fluid) THEN
                  ug( n2 ) = -ug( n1 )
                  us(n2,:) = -us(n1,:)
                END IF
                !
                IF(.NOT.BTEST(flag(ijmkp),0)) THEN
                !IF(flag(ijmkp) /= fluid) THEN
                  wg( n2 ) = -wg( n1 )
                  ws(n2,:) = -ws(n1,:)
                END IF
!
              CASE (free_io, zero_grad)
!  	    
                d0 = dy(j+1)
                d1 = dy(j)
                d2 = dy(j-1)
  
                !
                ! ... Normal velocities
                vg(n2) = vg(n1)
                vs(n2,:) = vs(n1,:)
                ! ... Tangential velocities
                !
                ug(n2)   = ug(n1)
                us(n2,:) = us(n1,:)
                wg(n2)   = wg(n1)
                ws(n2,:) = ws(n1,:)
                !
                ! ... Outlet scalars
                gnorm = gravy * 0.5D0 * (dy(j-1) + dy(j))
                IF (flag(n2) == zero_grad) gnorm = 0.D0
                CALL inoutflow(vg(n2),vs(n2,:),n2,n1,k,gnorm)
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
              ! ... Normal velocities
              !
              wg(n1) = 0.D0
              ws(n1,:) = 0.D0
              wg(n1) = 0.D0
              ws(n1,:) = 0.D0

              ! ... Tangential velocities
              !
              ug(n2)   = ug(n1)
              us(n2,:) = us(n1,:)
              IF (job_type == JOB_TYPE_3D) THEN
                vg(n2)   = vg(n1)
                vs(n2,:) = vs(n1,:)
              END IF
!
            CASE (noslip_wall)
!  
              ! ... Normal velocity
              wg(n1) = 0.D0
              ws(n1,:) = 0.D0
              wg(n1) = 0.D0
              ws(n1,:) = 0.D0

              ! ... Tangential velocities
              !
              IF (.NOT.BTEST(flag(ipjkp),0)) THEN
              !IF (flag(ipjkp) /= fluid) THEN
                ug(n2)   = -ug(n1)
                us(n2,is) = -us(n1,is)
              END IF

              IF (.NOT.BTEST(flag(ijpkp),0) .AND. job_type == JOB_TYPE_3D) THEN
              !IF (flag(ijpkp) /= fluid .AND. job_type == JOB_TYPE_3D) THEN
                vg(n2)   = -vg(n1)
                vs(n2,is) = -vs(n1,is)
              END IF
!
            CASE (free_io, zero_grad) 
!
              d0 = dz(k-1)
              d1 = dz(k)
              d2 = dz(k+1)

              ! ... Tangential velocities
              ug(n2)   = ug(n1)
              us(n2,:) = us(n1,:)
              IF (job_type == JOB_TYPE_3D) THEN
                vg(n2)   = vg(n1)
                vs(n2,:) = vs(n1,:)
              END IF

              ! ... Extrapolate the normal component of the velocities
              !
              wg(ijk)    = extrapolate(wg(ijkm), wg(ijkmm), d0, d1, ord)
              DO is = 1, nsolid
                ws(ijk,is)  = extrapolate(ws(ijkm,is), ws(ijkmm,is), d0, d1, ord)
              END DO
              wg(ijkp)   = extrapolate(wg(ijk), wg(ijkm), d1, d2, ord)
              ws(ijkp,:)   = extrapolate(ws(ijk,is), ws(ijkm,is), d1, d2, ord)
!
              ! ... Outlet scalars
              IF (flag(n2) == zero_grad) THEN
                CALL zerograd(n2,n1,k+1)
              ELSE IF (flag(n2) == free_io) THEN
                gnorm = gravz * 0.5D0 * (dz(k)+dz(k+1))
                CALL inoutflow(wg(n1),ws(n1,:),n2,n1,k+1,gnorm)
              END IF
!
            CASE DEFAULT
!
              CONTINUE

            END SELECT
!
! ***** Bottom boundary conditions ***** !
!
            n1 = ijk
            n2 = ijkm

            SELECT CASE ( flag( n2 ) )
!
            CASE (slip_wall)
!
              ! ... Normal velocities
              !
              wg( n2 ) = 0.D0
              ws(n2,:) = 0.D0

              ! ... Tangential velocities
              !
              IF(flag(ipjkm) == slip_wall) THEN 
                ug( n2 ) = ug( n1 )
                us(n2,:) = us(n1,:)
              END IF

              IF (flag(ijpkm) == slip_wall .AND. job_type == JOB_TYPE_3D) THEN
                vg( n2 ) = vg( n1 )
                vs(n2,:) = vs(n1,:)
              END IF
!
! ... Update sink term
!
              IF (isink > 0) THEN
                DO is = 1, nsolid
                   CALL sink_update(n1,is)
                END DO
              END IF
! ...
!
            CASE (noslip_wall)
!
              ! ... Normal velocities
              !
              wg(n2)   = 0.D0
              ws(n2,:) = 0.D0
!                
              IF(flag(ipjkm) == noslip_wall) THEN 
                ug(n2)   = -ug(n1)
                us(n2,:) = -us(n1,:)
              END IF

              IF(flag(ijpkm) == noslip_wall .AND. job_type == JOB_TYPE_3D) THEN
                vg(n2)   = -vg(n1)
                vs(n2,:) = -vs(n1,:)
              END IF
!
! ... Update sink term
!
              IF (isink > 0) THEN
                DO is = 1, nsolid
                   CALL sink_update(n1,is)
                END DO
              END IF
!
! ... The sink array must be updated also if the bottom cell is filled1 and filled2
!
            CASE (filled_cell_1, filled_cell_2)
              IF (isink > 0) THEN
                DO is = 1, nsolid
                   CALL sink_update(n1,is)
                END DO
              END IF
!
              CASE (free_io, zero_grad)
!  	    
                d0 = dz(k+1)
                d1 = dz(k)
                d2 = dz(k-1)
  
                !
                ! ... Normal velocities
                wg(n2) = wg(n1)
                ws(n2,:) = ws(n1,:)
                ! ... Tangential velocities
                !
                ug(n2)   = ug(n1)
                us(n2,:) = us(n1,:)
                vg(n2)   = vg(n1)
                vs(n2,:) = vs(n1,:)
                !
                ! ... Outlet scalars
                gnorm = gravz * 0.5D0 * (dz(k-1)+dz(k))
                IF (flag(n2) == zero_grad) gnorm = 0.D0
                CALL inoutflow(wg(n2),ws(n2,:),n2,n1,k-1,gnorm)
!
            CASE DEFAULT
!
              CONTINUE

            END SELECT
!
        END IF
      END DO mesh_loop
!
! ... Write immersed boundaries data
!
      IF (lpr > 2 .AND. immb >= 1) THEN
        !
        nfptx = SIZE(fptx)
        ALLOCATE(tfptx(nfptx))
        tfptx(:) = fptx(:)%vel(1)
        CALL parallel_sum_real(tfptx(:),nfptx)
        !
        IF (job_type == JOB_TYPE_3D) THEN
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
          IF (job_type == JOB_TYPE_3D) THEN
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
        IF (job_type == JOB_TYPE_3D) DEALLOCATE(tfpty)
        DEALLOCATE(tfptz)
      END IF
!
      RETURN
      END SUBROUTINE boundary
!-----------------------------------------------------------------------
      END MODULE boundary_conditions
!-----------------------------------------------------------------------
