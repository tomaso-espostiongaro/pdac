!-----------------------------------------------------------------------
      MODULE vent_conditions
!
! ... Identify vent cells on the topography
! ... Set vent conditions within vent cells
! ... Modify vent conditions in time
!
!-----------------------------------------------------------------------
      USE dimensions, ONLY: max_nsolid, max_ngas, nsolid, ngas
      USE io_files, ONLY: errorunit, ventunit, ventfile
      IMPLICIT NONE
      PUBLIC
      !
      ! ... flags
      !
      INTEGER :: ivent, irand, iali, ipro
      !
      CHARACTER(LEN=80) :: rad_file
      !
      REAL*8 :: xvent, yvent, vent_radius, base_radius, crater_radius
      REAL*8 :: wrat
      REAL*8 :: u_gas, v_gas, w_gas, p_gas, t_gas
      REAL*8 :: u_solid(max_nsolid), v_solid(max_nsolid), w_solid(max_nsolid), &
                ep_solid(max_nsolid), t_solid(max_nsolid)
      REAL*8 :: vent_ygc(max_ngas)
!      
      REAL*8, ALLOCATABLE, DIMENSION(:) :: rad
      REAL*8, ALLOCATABLE, DIMENSION(:) :: ug_rad, wg_rad, &
                                           p_rad, tg_rad
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: ygc_rad
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: us_rad, ws_rad, &
                                             ep_rad, ts_rad
      TYPE icvent_cell
        INTEGER :: imesh
        INTEGER :: i
        INTEGER :: j
        INTEGER :: k
        REAL*8  :: frac
        REAL*8  :: fact
      END TYPE icvent_cell

      TYPE(icvent_cell), ALLOCATABLE :: vcell(:)

      INTEGER :: nvt, iiv, jjv, kkv
      INTEGER :: seed
!
      PRIVATE :: icvent_cell, vcell, nvt, seed
      PRIVATE :: ug_rad, wg_rad, p_rad, tg_rad, rad, ygc_rad,  &
                 us_rad, ws_rad, ep_rad, ts_rad
      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
! ... This routine locate the vent on a georeferenced mesh and set the
! ... cell flags onto the neighbours of the inlet
!
      SUBROUTINE locate_vent

      USE control_flags, ONLY: job_type, lpr
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE grid, ONLY: x, y, z, iv, jv
      USE parallel, ONLY: mpime, root

      IMPLICIT NONE
      INTEGER :: n
!
! ... If the vent coordinates are not assigned,
! ... the mesh center coordinates are assumed
!
      IF (xvent == -99999.D0) xvent = x(iv)
      IF (yvent == -99999.D0) yvent = y(jv)
!
      IF (job_type == JOB_TYPE_2D ) THEN
              CALL locate_vent_2D
      ELSE IF (job_type == JOB_TYPE_3D) THEN
              CALL locate_vent_3D
      END IF
!
! ... Print out the vent coordinates on the standard output
! ... and the vent conditions
!
      IF( lpr > 0 .AND. mpime == root ) THEN
        OPEN(UNIT=ventunit,FILE=ventfile,STATUS='UNKNOWN')
        WRITE(ventunit,100) iiv, jjv, kkv
        WRITE(ventunit,200) x(iiv), y(jjv), z(kkv)
        WRITE(ventunit,300) vent_radius, base_radius
        WRITE(ventunit,*) 
        WRITE(ventunit,*) 'Number of vent cells: ', nvt
        WRITE(ventunit,*) 'Vent cells report: '
        DO n = 1, nvt
          WRITE(ventunit,400) n, vcell(n)
        END DO
      END IF
!
 100    FORMAT(1X,'vent center: ',3I5)
 200    FORMAT(1X,'vent center coordinates: ',3(F12.2))
 300    FORMAT(1X,'vent radius and base radius: ',2(F12.2))
 400    FORMAT(I4,4(I5),2(F8.4))
!
      RETURN
      END SUBROUTINE locate_vent
!-----------------------------------------------------------------------
! ... This routine locate the vent on a georeferenced mesh and set the
! ... cell flags onto the neighbours of the inlet
!
      SUBROUTINE locate_vent_2D

      USE atmospheric_conditions, ONLY: p_atm
      USE control_flags, ONLY: job_type, lpr
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nx, ny, nz
      USE grid, ONLY: x, z, fl, xb, zb, dz, dx, itc
      USE grid, ONLY: fluid, vent_cell, noslip_wall, slip_wall
      USE volcano_topography, ONLY: itp, ord
      USE volcano_topography, ONLY: flatten_dem_vent_2D
      USE volcano_topography, ONLY: rim_quota
      USE parallel, ONLY: mpime, root
!
      IMPLICIT NONE
      
      INTEGER :: i, j, k
      INTEGER :: ijk, nv
      INTEGER :: iwest, ieast, jnorth, jsouth
      INTEGER :: quota, dk
      REAL*8 :: soglia
!
! ... Vent coordinates '(xvent, yvent)' are taken from the input.
! ... Notice that the vent center coincides with the cell left edge
! ... (different from 3D)
!
      DO i = 1, nx
        IF (x(i) <= xvent) iiv = i
      END DO
      IF (iiv == 1 .OR. itc == 1) iiv = 2
      xvent = xb(iiv-1)
!
! ... Define the 'quota' of the volcanic vent
! ... (considering the topography). If the topography
! ... has been red from a file, 'kkv' has already been set
! ... in the 'volcano_topography' module.
!
      IF( itp < 1 ) THEN
        DO k= 1, nz
          ijk = iiv + (k-1) * nx 
          IF (fl(ijk) == slip_wall .OR. fl(ijk) == noslip_wall) kkv = k
        END DO
        quota = kkv
      ELSE
        quota = ord(iiv)
        !
        ! ... Modify the topography at the base or at the top
        ! ... of the crater to flatten the profile around the vent. 
        ! ... Reset the cell flags and the 'dist' array.
        !
        CALL flatten_dem_vent_2D(xvent,vent_radius,quota)
      END IF
!
! ... Define the rectangle enclosing the vent.
! ... 'nvt' is the number of vent cells
!
      iwest = 2
      DO i = 2, nx
        IF (xb(i-1) <= (xvent-vent_radius)) iwest = i
        IF (xb(i-1) < (xvent+vent_radius)) ieast = i
      END DO
      nvt = (ieast-iwest+1)
!
! ... allocate and initialize the cells enclosing the vent
!
      ALLOCATE(vcell(nvt))
      vcell(:)%frac = 1.D0
      vcell(:)%fact = 1.D0
!      
! ... Loop over the cells enclosing the vent
!
      nv = 0
      DO i = iwest, ieast
        nv = nv + 1

        ! ... Below the vent quota, set the cell flag as noslip walls
        !
        DO k = 1, quota - 1
          ijk = i + (k-1) * nx
          fl(ijk) = noslip_wall
        END DO

        ! ... At the vent quota, vent cell flags are set to vent_cell
        !
        k = quota
        ijk = i + (k-1) * nx
        
        vcell(nv)%imesh = ijk
        vcell(nv)%i = i
        vcell(nv)%j = 0
        vcell(nv)%k = k
        vcell(nv)%frac = 1.D0
        fl(ijk) = vent_cell
        !
        IF (nv == nvt) THEN
                IF (itc == 1) THEN
                        vcell(nvt)%frac = 1.D0 - (xb(ieast) - xb(iwest-1) - vent_radius)/ dx(ieast)
                ELSE IF (itc == 0) THEN
                        vcell(nvt)%frac = 1.D0 - (xb(ieast) - xb(iwest-1) - 2.D0*vent_radius)/ dx(ieast)
                END IF
        END IF
        
        ! ... Above the vent quota, cell flags are set to '1'
        ! ... (fluid cells)
        !
        DO k = quota+1, nz-1
          ijk = i + (k-1) * nx
          fl(ijk) = fluid
        END DO
                  
      END DO
!
      RETURN
!
      END SUBROUTINE locate_vent_2D
!-----------------------------------------------------------------------
! ... This routine locate the vent on a georeferenced mesh and set the
! ... cell flags onto the neighbours of the inlet
!
      SUBROUTINE locate_vent_3D

      USE atmospheric_conditions, ONLY: p_atm
      USE control_flags, ONLY: job_type, lpr
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nx, ny, nz
      USE grid, ONLY: x, y, z, fl, xb, yb, zb, dz
      USE grid, ONLY: fluid, vent_cell, noslip_wall, slip_wall
      USE volcano_topography, ONLY: itp, iavv, ord2d
      USE volcano_topography, ONLY: nocrater, flatten_dem_vent
      USE volcano_topography, ONLY: rim_quota
      USE parallel, ONLY: mpime, root

      IMPLICIT NONE
      
      INTEGER :: i, j, k
      INTEGER :: ijk, nv
      INTEGER :: iwest, ieast, jnorth, jsouth
      INTEGER :: quota, dk
      REAL*8 :: soglia
!
      ! ... The crater base radius must include, at least, 
      ! ... the rectancle containing the vent
      !
      IF( base_radius < 1.6D0*vent_radius ) THEN
              base_radius = 1.6D0*vent_radius
              IF (mpime == root) &
                WRITE(errorunit,*) 'WARNING! control the crater base!'
      END IF
!
! ... The coordinates '(xvent, yvent)' are taken from the input.
! ... Notice that the vent center coincides with the cell center 
! ... (different from 2D)
!
      DO i = 1, nx
        IF (x(i) <= xvent) iiv = i
      END DO
      DO j = 1, ny
        IF (y(j) <= yvent) jjv = j
      END DO
      xvent = x(iiv)
      yvent = y(jjv)
!
! ... Define the 'quota' of the volcanic vent
! ... (considering the topography). If the topography
! ... has been red from a file, 'kkv' has already been set
! ... in the 'volcano_topography' module.
!
      IF( itp < 1 ) THEN

        DO k= 1, nz
          ijk = iiv + (jjv-1) * nx + (k-1) * nx * ny
          IF (fl(ijk) == slip_wall .OR. fl(ijk) == noslip_wall) kkv = k
        END DO
        quota = kkv

      ELSE
!
        kkv = ord2d(iiv,jjv)
        !
        ! ... Reset the vent quota if the topography has to be modified
        !
        IF (nocrater) THEN
          DO k = 1, nz
            IF (zb(k) <= rim_quota) kkv = k
          END DO
        END IF
        quota = kkv
!
        ! ... Modify the topography at the base or at the top
        ! ... of the crater to flatten the profile around the vent. 
        ! ... Reset the cell flags and the 'dist' array.
        !
        CALL flatten_dem_vent(xvent,yvent,base_radius,crater_radius,quota)
!
      END IF
!
! ... Define the rectangle enclosing the vent.
! ... 'nvt' is the number of vent cells
!
      iwest = 1
      jsouth = 1
      DO i = 2, nx
        IF (xb(i-1) <= (xvent-vent_radius)) iwest = i
        IF (xb(i-1) < (xvent+vent_radius)) ieast = i
      END DO
      DO j = 2, ny
        IF (yb(j-1) <= (yvent-vent_radius)) jsouth = j
        IF (yb(j-1) < (yvent+vent_radius)) jnorth = j
      END DO
      nvt = (ieast-iwest+1)*(jnorth-jsouth+1)
!
! ... allocate and initialize the cells enclosing the vent
!
      ALLOCATE(vcell(nvt))
      vcell(:)%frac = 1.D0
      vcell(:)%fact = 1.D0
!
! ... Set the threshold above which a cells is included
! ... into the vent
!
      IF (irand == 0 .AND. iali == 0) THEN
        soglia = 0.5D0
      ELSE
        soglia = 0.D0
      END IF
!      
! ... Loop over the cells enclosing the vent
!
      nv = 0
      DO j = jsouth, jnorth
        DO i = iwest, ieast
          nv = nv + 1

          ! ... Below the vent quota, set the cell flag as noslip walls
          !
          DO k = 1, quota - 1
            ijk = i + (j-1) * nx + (k-1) * nx * ny
            fl(ijk) = noslip_wall
          END DO

          ! ... At the vent quota, vent cell flags are set to vent_cell
          ! ... The fraction of the cells occupied by the vent
          ! ... is computed.
          !
          k = quota
          ijk = i + (j-1) * nx + (k-1) * nx * ny
          
          vcell(nv)%imesh = ijk
          vcell(nv)%i = i
          vcell(nv)%j = j
          vcell(nv)%k = k
          vcell(nv)%frac = cell_fraction(i,j)
          
          IF (vcell(nv)%frac > soglia) THEN
            fl(ijk) = vent_cell
          ELSE
            fl(ijk) = noslip_wall
          END IF
          
          ! ... Above the vent quota, cell flags are set to '1'
          ! ... (fluid cells)
          !
          DO k = quota+1, nz-1
            ijk = i + (j-1) * nx + (k-1) * nx * ny
            fl(ijk) = fluid
          END DO
                    
        END DO
      END DO
!
      RETURN
      END SUBROUTINE locate_vent_3D
!-----------------------------------------------------------------------
      SUBROUTINE set_ventc
!
! ... Compute the steady inlet conditions for a circular vent
! ... An aribtrary radial profile can be assigned, as a function
! ... of the averaged vertical velocity
!
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nsolid, ngas
      USE domain_mapping, ONLY: ncint, meshinds
      USE environment, ONLY: cpclock
      USE eos_gas, ONLY: ygc
      USE gas_constants, ONLY: gas_type
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_temperature, ONLY: tg, ts
      USE gas_solid_velocity, ONLY: ug, wg, vg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE grid, ONLY: flag, x, y, vent_cell
      USE io_files, ONLY: testunit
      USE parallel, ONLY: mpime, root
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: ep, p
      USE array_filters, ONLY: interp
      IMPLICIT NONE

      REAL*8 :: ygcsum, mfr
      INTEGER :: ijk, imesh, i,j,k, is, ig, n
      REAL*8 :: alpha, beta, ra, dex, dey, angle, angle4
      REAL*8 :: fact_r
      REAL*8 :: distance
      REAL*8 :: rseed
!
! ... Compatibility
      IF (ipro >= 1) THEN
            CALL read_radial_profile
            iali = 0
            irand = 0
            wrat = 1.D0
      END IF
      IF (irand == 1) THEN
              iali = 0
              wrat = 1.D0
      END IF
! 
      DO ijk = 1, ncint      
        IF(flag(ijk) == vent_cell) THEN
          CALL meshinds(ijk,imesh,i,j,k)
!
          ! ... Determine the fraction of the cell
          ! ... not occupied by the topography.
          !
          DO n = 1, nvt
            IF (vcell(n)%imesh == imesh) THEN
              alpha = vcell(n)%frac
            END IF
          END DO
!          
          IF (ipro >= 1) THEN
            !
            IF (job_type == JOB_TYPE_2D ) THEN
              distance = (x(i) - xvent)
              angle = 0.D0
            ELSE IF (job_type == JOB_TYPE_3D) THEN
              distance = DSQRT( (x(i)-xvent)**2 + (y(j)-yvent)**2 )
              angle    = ATAN2( (y(j)-yvent), (x(i)-xvent) )
            END IF
            !
            CALL interp(rad, ug_rad, distance, ug(ijk))
            IF (job_type == JOB_TYPE_3D) vg(ijk) = ug(ijk) * SIN(angle)
            ug(ijk) = ug(ijk) * COS(angle)
            CALL interp(rad, wg_rad, distance, wg(ijk))
            !
            CALL interp(rad, tg_rad, distance, tg(ijk))
            CALL interp(rad, p_rad,  distance, p(ijk))
            !
            DO ig = 1, ngas
              CALL interp(rad, ygc_rad(:,ig), distance, ygc(ijk,ig))
            END DO
            !
            ep(ijk) = 1.D0
            DO is = 1,nsolid
              CALL interp(rad, us_rad(:,is),  distance, us(ijk,is))
              IF (job_type == JOB_TYPE_3D) vs(ijk,is) = us(ijk,is) * SIN(angle)
              us(ijk,is) = us(ijk,is) * COS(angle)
              CALL interp(rad, ws_rad(:,is),  distance, ws(ijk,is))
              !
              CALL interp(rad, ts_rad(:,is),  distance, ts(ijk,is))
              CALL interp(rad, ep_rad(:,is), distance, rlk(ijk,is))
              ep(ijk) = ep(ijk) - rlk(ijk,is)
              rlk(ijk,is) = rlk(ijk,is) * rl(is)
            END DO
            !
          ELSE
            ! ... Set uniform initial conditions on
            ! ... all cells enclosing the vent,
            ! ... as specified in the input namelist 'inlet'
            !
            ug(ijk) = u_gas 
            IF (job_type == JOB_TYPE_3D) vg(ijk) = v_gas 
            wg(ijk) = w_gas
            !
            tg(ijk) = t_gas
            ep(ijk) = 1.D0 - SUM(ep_solid(1:nsolid))
            p(ijk)  = p_gas
            !
            DO ig = 1, ngas
              ygc(ijk,ig) = vent_ygc(gas_type(ig))
            END DO
            !
            DO is = 1,nsolid
              us(ijk,is)  = u_solid(is) 
              IF (job_type == JOB_TYPE_3D) vs(ijk,is)  = v_solid(is) 
              ws(ijk,is) = w_solid(is)
              !
              ts(ijk,is)  = t_solid(is)
              rlk(ijk,is) = ep_solid(is)*rl(is)
            END DO
          END IF
!          
          ! ... check gas components closure relation
          !
          ygcsum = SUM(ygc(ijk,:))
          IF ( ygcsum /= 1.D0 ) THEN
            ygc(ijk,ngas) = 1.D0 - SUM( ygc(ijk,1:ngas-1) )
          END IF
!          
          ! ... Mixture density/velocity is corrected in those cells
          ! ... partially filled by the topography in order
          ! ... to respect the mass flux
          !
          IF (iali == 1) CALL density_antialias(ijk,k,alpha)
          IF (iali == 2) CALL density_antialias_2(ijk,k,alpha)
          IF (iali == 3) CALL velocity_antialias(ijk,alpha)

          ! ... 'wrat' is the ratio between the maximum
          ! ... vertical velocity and the averaged velocity
          ! ... If 'wrat' is greater than 1, the inlet profile
          ! ... decreases to 0 towards the vent rim as a power law,
          ! ... and the surface average of the velocity equals the 
          ! ... input velocity...
          ! 
          IF (wrat > 1.D0 .AND. iali /= 3) THEN
            ! ... linear average 
            !beta = 1.D0 / (wrat - 1.D0)
            ! ... surface average 
            beta = 2.D0 / (wrat - 1.D0)
            dey = y(j)-yvent
            dex = x(i)-xvent
            ra = DSQRT(dex**2 + dey**2)
            ra = MIN(ra / vent_radius, 1.D0)
            fact_r = wrat * (1.D0 - ra ** beta)
            !
            ! ... Angular modulation of the velocity profile
            !
            ! IF (dey /= 0.D0 ) THEN
            !   angle = ATAN2(dey,dex)
            ! ELSE IF (dex >= 0.D0) THEN
            !   angle = 0.D0
            ! ELSE IF (dex < 0.D0) THEN
            !   angle = 4.D0 * ATAN(1.D0)
            ! END IF
            ! angle4 = angle*4.D0
            ! fact_r = fact_r + 0.1D0 * COS(angle4)
            !
            CALL correct_velocity_profile(ijk,fact_r)
          END IF

          ! ... determine the initial random seed
          !
          !CALL MP_WALLTIME(rseed,mpime)
          rseed = cpclock()
          seed = INT(rseed)
          WRITE(testunit,*) 'seed=', seed
          !IF (mpime == root) seed = INT(cpclock())
          !CALL bcast_integer(seed,1,root)
        
        END IF
      END DO

      RETURN
      END SUBROUTINE set_ventc
!-----------------------------------------------------------------------
      SUBROUTINE read_radial_profile
      USE parallel, ONLY: mpime, root
      USE io_files, ONLY: tempunit
      IMPLICIT NONE
      INTEGER :: raddim, is, n, ig
!
      IF (mpime == root) THEN
              OPEN(tempunit,FILE=rad_file,STATUS='OLD')
              READ(tempunit,*) raddim
      END IF
!
      CALL bcast_integer(raddim,1,root)
!
      ALLOCATE (ug_rad(raddim),wg_rad(raddim), &
          tg_rad(raddim), p_rad(raddim), rad(raddim))
      ALLOCATE (ygc_rad(raddim,ngas))
      ALLOCATE (us_rad(raddim,nsolid), &
          ws_rad(raddim,nsolid), ts_rad(raddim,nsolid), ep_rad(raddim,nsolid))
!
      IF (mpime == root) THEN
              READ(tempunit,*) (rad(n), n=1, raddim)
              READ(tempunit,*) (ug_rad(n), n=1, raddim)
              READ(tempunit,*) (wg_rad(n), n=1, raddim)
              READ(tempunit,*) (tg_rad(n), n=1, raddim)
              READ(tempunit,*) (p_rad(n), n=1, raddim)
              DO ig = 1, ngas
                READ(tempunit,*) (ygc_rad(n,ig), n=1, raddim)
              END DO
              DO is = 1, nsolid
                READ(tempunit,*) (us_rad(n,is), n=1, raddim)
                READ(tempunit,*) (ws_rad(n,is), n=1, raddim)
                READ(tempunit,*) (ts_rad(n,is), n=1, raddim)
                READ(tempunit,*) (ep_rad(n,is), n=1, raddim)
              END DO
              CLOSE(tempunit)
      END IF
!
      CALL bcast_real(rad,raddim,root)
      CALL bcast_real(ug_rad,raddim,root)
      CALL bcast_real(wg_rad,raddim,root)
      CALL bcast_real(tg_rad,raddim,root)
      CALL bcast_real(p_rad,raddim,root)
      CALL bcast_real(ygc_rad,raddim*ngas,root)
      CALL bcast_real(us_rad,raddim*nsolid,root)
      CALL bcast_real(ws_rad,raddim*nsolid,root)
      CALL bcast_real(ts_rad,raddim*nsolid,root)
      CALL bcast_real(ep_rad,raddim*nsolid,root)
!
      RETURN
!
      END SUBROUTINE read_radial_profile
!-----------------------------------------------------------------------
      SUBROUTINE density_antialias(ijk,k,alpha)
!
! ... Correct the density in the inlet cells partially filled by the
! ... topography ("antialiasing")
!
      USE atmospheric_conditions, ONLY: p_atm, t_atm
      USE dimensions, ONLY: nsolid, ngas
      USE eos_gas, ONLY: xgc, ygc, mole
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_temperature, ONLY: tg
      USE gas_constants, ONLY: gmw, gas_type
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: ep, p
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ijk, k
      REAL*8, INTENT(IN) :: alpha
      REAL*8 :: ep0, mg, rrhoatm
      INTEGER :: is, ig
      
      ! ... 'alpha' is the cell fraction.
      ! ... The solid bulk density is reduced by the factor 'alpha'
      !
      DO is = 1,nsolid
        rlk(ijk,is) = ep_solid(is)*rl(is) * alpha
      END DO
      
      ! ... The density contrast is also reduced by the factor 'alpha'
      ! ... by correcting the void fraction and the pressure.
      !
      ep0     = 1.D0 - SUM(ep_solid(1:nsolid))
      ep(ijk) = 1.D0 - alpha * SUM(ep_solid(1:nsolid))

      CALL mole( xgc(ijk,:), ygc(ijk,:) )

      mg = 0.D0
      DO ig = 1, ngas
         mg = mg + xgc(ijk,ig) * gmw(gas_type(ig))
      END DO

      rrhoatm = p_atm(k) / t_atm(k) * gmw(6)
      p(ijk) = p_gas * alpha * ep0 / ep(ijk) + &
               rrhoatm / mg * tg(ijk) * (1.D0-alpha) / ep(ijk)

      RETURN
      END SUBROUTINE density_antialias
!-----------------------------------------------------------------------
      SUBROUTINE density_antialias_2(ijk,k,alpha)
!
! ... Correct the density in the inlet cells partially filled by the
! ... topography ("antialiasing")
!
      USE eos_gas, ONLY: xgc, mole, thermal_eosg
      USE dimensions, ONLY: nsolid
      USE gas_constants, ONLY: gmw, gas_type
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_temperature, ONLY: tg
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: ep, p
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ijk, k
      REAL*8, INTENT(IN) :: alpha
      REAL*8 :: ep0, rog0, rlk0, beta
      INTEGER :: is, ig
      
      CALL thermal_eosg(rog0, tg(ijk), p(ijk), xgc(ijk,:))
      ep0  = 1.D0 - SUM(ep_solid(1:nsolid))
      rlk0 = SUM(rlk(ijk,:))
      beta = rlk0 / (rlk0 + rog0*ep0)
      !
      ! ... 'alpha' is the cell surface fraction.
      ! ... 'beta' is the total solid mass fraction
      ! ... The solid bulk density is reduced by the factor 'alpha/beta'
      !
      DO is = 1,nsolid
        rlk(ijk,is) = ep_solid(is)*rl(is) * alpha / beta
      END DO
      ep(ijk) = 1.D0 - alpha/beta * SUM(ep_solid(1:nsolid))
      
      RETURN
      END SUBROUTINE density_antialias_2
!-----------------------------------------------------------------------
      SUBROUTINE velocity_antialias(ijk,alpha)
!
! ... Correct the velocity in the inlet cells partially filled by the
! ... topography ("antialiasing")
!
      USE dimensions, ONLY: nsolid
      USE gas_solid_velocity, ONLY: wg, ws
      IMPLICIT NONE

      REAL*8, INTENT(IN) :: alpha
      INTEGER, INTENT(IN) :: ijk
      INTEGER :: is

      wg(ijk) = w_gas * alpha
      
      DO is = 1,nsolid
        ws(ijk,is) = w_solid(is) * alpha
      END DO

      RETURN
      END SUBROUTINE velocity_antialias
!-----------------------------------------------------------------------
      SUBROUTINE correct_velocity_profile(ijk,factor)
!
! ... Compute the steady inlet conditions for a circular vent
! ... An aribtrary radial profile can be assigned, as a function
! ... of the averaged vertical velocity
!
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nsolid
      USE gas_solid_velocity, ONLY: ug, us, vg, vs, wg, ws
      USE pressure_epsilon, ONLY: p
      IMPLICIT NONE

      REAL*8, INTENT(IN) :: factor
      INTEGER, INTENT(IN) :: ijk
      INTEGER :: is

      ug(ijk) = u_gas * factor
      IF (job_type == JOB_TYPE_3D) vg(ijk) = v_gas * factor
      wg(ijk) = w_gas * factor
      
      DO is = 1,nsolid
        us(ijk,is) = u_solid(is) * factor
        IF (job_type == JOB_TYPE_3D) vs(ijk,is) = v_solid(is) * factor
        ws(ijk,is) = w_solid(is) * factor
      END DO

      RETURN
      END SUBROUTINE correct_velocity_profile
!-----------------------------------------------------------------------
      SUBROUTINE update_vent_cell(ijk,imesh,sweep)
!
      USE dimensions, ONLY: nsolid
      USE gas_solid_velocity, ONLY: wg, ws
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ijk, imesh, sweep
      INTEGER :: is, n
      REAL*8 :: switch
      
      ! ... The inlet cell is on/off accordingly to the
      ! ... value 1/0 of 'switch'
      ! 
      IF (irand == 1) THEN
        DO n = 1, nvt
          IF (vcell(n)%imesh == imesh) THEN
            switch = vcell(n)%fact
          END IF
        END DO
      ELSE
        switch = 1.D0
      END IF

      wg(ijk) = w_gas * switch
      DO is = 1,nsolid
        ws(ijk,is) = w_solid(is) * switch
      END DO

      RETURN
      END SUBROUTINE update_vent_cell
!-----------------------------------------------------------------------
      SUBROUTINE update_inlet_cell(ijk)
!
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nsolid
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      IMPLICIT NONE

      REAL*8 :: ran0
      EXTERNAL :: ran0
      INTEGER, INTENT(IN) :: ijk
      INTEGER :: is, n
      
      ! ... The velocities in a cell are perturbed accordingly
      ! ... to a random function 'ran0' with a maximum intensity
      ! ... of 1%
      !
      ug(ijk) = ug(ijk) * (1.D0 + 0.01D0 * (ran0(seed) - 0.5D0) )
      IF (job_type == JOB_TYPE_3D)  &
        vg(ijk) = vg(ijk) * (1.D0 + 0.01D0 * (ran0(seed) - 0.5D0) )
      wg(ijk) = wg(ijk) * (1.D0 + 0.01D0 * (ran0(seed) -0.5D0) )
      
      DO is = 1,nsolid
        us(ijk,is) = us(ijk,is) * (1.D0 + 0.01D0 * (ran0(seed) -0.5D0) )
        IF (job_type == JOB_TYPE_3D) &
          vs(ijk,is) = vs(ijk,is) * (1.D0 + 0.01D0 * (ran0(seed) -0.5D0) )
        ws(ijk,is) = ws(ijk,is) * (1.D0 + 0.01D0 * (ran0(seed) -0.5D0) )
      END DO

      RETURN
      END SUBROUTINE update_inlet_cell
!-----------------------------------------------------------------------
      SUBROUTINE grow_vent_cell(ijk,imesh,sweep)
!
      USE dimensions, ONLY: nsolid
      USE gas_solid_velocity, ONLY: wg, ws
      USE time_parameters, ONLY: ift
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ijk, imesh, sweep
      INTEGER :: is, n
      REAL*8 :: switch, growth_factor
      
      IF (wrat > 1.D0) RETURN
      
      ! ... The inlet profile grows from 0 to its final
      ! ... value as the function of time
      !
      SELECT CASE (ift)
        CASE (1)
          growth_factor = ft1(sweep)
        CASE (2)
          growth_factor = ft2(sweep)
        CASE (3)
          growth_factor = ft3(sweep)
        CASE (4)
          growth_factor = ft4(sweep)
        CASE (5)
          growth_factor = ft5(sweep)
        CASE DEFAULT
          growth_factor = 1.D0
      END SELECT
!
      wg(ijk) = w_gas * growth_factor
      DO is = 1,nsolid
        ws(ijk,is) = w_solid(is) * growth_factor
      END DO

      RETURN
      END SUBROUTINE grow_vent_cell
!-----------------------------------------------------------------------
      SUBROUTINE grow_inlet_cell(ijk,imesh,sweep)
!
      USE atmospheric_conditions, ONLY: p_ground
      USE dimensions, ONLY: nsolid
      USE gas_solid_velocity, ONLY: wg, ws
      USE time_parameters, ONLY: ift
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ijk, imesh, sweep
      INTEGER :: is, n
      REAL*8 :: growth_factor
      REAL*8 :: wgas_init, wpart_init(max_nsolid)
      
      IF (sweep == 1) THEN
              wgas_init = wg(ijk) 
              wpart_init(1:nsolid) = ws(ijk,1:nsolid)
      END IF
      !
      ! ... The inlet profile grows from 0 to its final
      ! ... value as the function of time
      !
      SELECT CASE (ift)
        CASE (1)
          growth_factor = ft1(sweep)
        CASE (2)
          growth_factor = ft2(sweep)
        CASE (3)
          growth_factor = ft3(sweep)
        CASE (4)
          growth_factor = ft4(sweep)
        CASE (5)
          growth_factor = ft5(sweep)
        CASE DEFAULT
          growth_factor = 1.D0
      END SELECT
      !
      wg(ijk) = wgas_init * growth_factor
      DO is = 1,nsolid
        ws(ijk,is) = wpart_init(is) * growth_factor
      END DO

      RETURN
      END SUBROUTINE grow_inlet_cell
!-----------------------------------------------------------------------
      SUBROUTINE random_switch(sweep)
!
! ... Randomly switch vent cells on/off, with a probability
! ... equal to the cell fraction
!
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: sweep
      REAL*8 :: rnv
      REAL*8 :: ran0
      EXTERNAL :: ran0
      INTEGER :: n

      IF (job_type == JOB_TYPE_2D ) RETURN

      DO n = 1, nvt
        rnv = ran0(seed)
        IF (rnv <= vcell(n)%frac) THEN
          vcell(n)%fact = 1.D0
        ELSE
          vcell(n)%fact = 0.D0
        END IF
      END DO

      RETURN
      END SUBROUTINE random_switch
!-----------------------------------------------------------------------
      REAL*8 FUNCTION ft1(n)
      USE time_parameters, ONLY: tau1, tau2, dt
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL*8 :: t, tr1, tr2

      t = n * dt
      tr1 = 0.D0
      tr2 = 0.D0
      ft1 = 1.D0
!
! ... The function 'ft1' grows linearly from 0 to 1 in a time 'tau1'
! ... and decreases linearly from 1 to 0 in a time 'tau2-tau1'
!
      IF (tau1 > 0.D0) tr1 = t / tau1
      IF (tau2 > tau1) tr2 = (tau2 - t)/(tau2 - tau1)

      IF (tr1 > 0.D0 .AND. tr2 > 0.D0) THEN
              ft1 = MAX(MIN(tr1,tr2),0.D0)
      END IF

      RETURN
      END FUNCTION ft1
!-----------------------------------------------------------------------
      REAL*8 FUNCTION ft2(n)
      USE time_parameters, ONLY: tau1, tau2, dt
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n        ! Dichiara l'argomento in ingresso
      REAL*8 :: t
! ... actual time
      t = n * dt
!    
      ft2 = DEXP(-(t-tau1)**2/(2*tau2**2))
!
      RETURN
      END FUNCTION ft2
!----------------------------------------------------------------------
      REAL*8 FUNCTION ft3(n)
      USE time_parameters, ONLY: tau1, tau2, dt
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n        ! Dichiara l'argomento in ingresso
      REAL*8 :: t
! ... actual time
      t = n * dt
    
! ... The function 'ft3' grows linearly from 0 to 1 in a time 'tau1'
! ... and decreases to 0 after time 'tau2'
!
      IF (t > 0.D0) THEN
        IF (t <= tau1) THEN
          ft3 = t / tau1   
        ELSE IF (t > tau1 .AND. t < 2*tau1 ) THEN
          ft3= -t/tau1 + 2.D0
        ELSE
          ft3 = 0.D0
        END IF
      END IF

      RETURN
      END FUNCTION ft3
!-----------------------------------------------------------------------
      REAL*8 FUNCTION ft4(n)
      USE time_parameters, ONLY: tau1, tau2, dt
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n        ! Dichiara l'argomento in ingresso
      REAL*8 :: t
      REAL*8 :: pi
!
      pi = 4.D0 * DATAN(1.D0)
!
! ... actual time
      t = n * dt
!    
! ... The function 'ft4' oscillates periodically from 0 to 1 with period 'tau1'
!
      IF (t > 0.D0 .AND. t <= tau2) THEN
        ft4 = (DSIN(t/tau1*pi))**2
      ELSE IF (t > tau2) THEN
        ft4 = 0.D0
      END IF

      RETURN
      END FUNCTION ft4
!-----------------------------------------------------------------------
      REAL*8 FUNCTION ft5(n)
      USE time_parameters, ONLY: tau1, tau2, dt
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL*8 :: t, tr1, tr2

      t = n * dt
      tr1 = 0.D0
      tr2 = 0.D0
!
! ... The function 'ft5' grows smoothly from 0 to 1 in a time 'tau1'
! ... and it is set to 0 after time 'tau2'
!
      IF (tau1 > 0.D0) tr1 = t / tau1
      IF (tau2 > 0.D0) tr2 = t / tau2

      IF (tr1 > 0.D0 .AND. tr2 > 0.D0) THEN
        IF (tr1 <= 1.D0) THEN
          ft5 = 3.D0 * (tr1)**2 - 2.D0 * (tr1)**3     
        ELSE IF (tr2 > 1.D0 ) THEN
          ft5 = 0.D0
        ELSE
          ft5 = 1.D0
        END IF
      END IF

      RETURN
      END FUNCTION ft5
!-----------------------------------------------------------------------
      REAL*8 FUNCTION cell_fraction(i,j)
!
! ... Assign a weight to partially filled vent cells:
! ... the weight is proportional to the area of the intersection of
! ... the cell and the vent 
!
      USE grid, ONLY: xb, yb, zb
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i, j
      
      TYPE polypoint
        REAL*8 :: x
        REAL*8 :: y
        INTEGER :: side             !(0=corner, 1=W, 2=N, 3=E, 4=S)
      END TYPE polypoint
      TYPE(polypoint) :: poly_p(16) !intersection points of bound.cell and vent 
                                    !(from west-south corner, clowckwise)
      INTEGER :: n,l
      
      REAL*8 :: delta_sol,x_sol1,x_sol2,y_sol1,y_sol2
      REAL*8 :: area_p
      
      n=0
      
      IF ( (xb(i-1)-xvent)**2 + (yb(j-1)-yvent)**2 <= vent_radius**2 ) THEN
         n=n+1
         poly_p(n)%x=xb(i-1)
         poly_p(n)%y=yb(j-1)
         poly_p(n)%side=0
      ENDIF
      
      delta_sol=vent_radius**2-(xb(i-1)-xvent)**2
      
      IF ( delta_sol > 0 ) THEN
         y_sol1= yvent-SQRT (delta_sol)
         y_sol2= yvent+SQRT (delta_sol)
         IF ( (yb(j-1) < y_sol1).AND. (y_sol1 < yb(j)) ) THEN
            n=n+1
            poly_p(n)%x=xb(i-1)
            poly_p(n)%y=y_sol1
            poly_p(n)%side=1
         ENDIF 
         
         IF ( (yb(j-1) < y_sol2) .AND. (y_sol2 < yb(j)) ) THEN
            n=n+1
            poly_p(n)%x=xb(i-1)
            poly_p(n)%y=y_sol2
            poly_p(n)%side=1
         ENDIF 
      ENDIF
      
      IF ( (xb(i-1)-xvent)**2 + (yb(j)-yvent)**2 <= vent_radius**2 ) THEN
         n=n+1
         poly_p(n)%x=xb(i-1)
         poly_p(n)%y=yb(j)
         poly_p(n)%side=0
      ENDIF
      
      delta_sol=vent_radius**2-(yb(j)-yvent)**2
      
      IF ( delta_sol > 0 ) THEN
         x_sol1= xvent-SQRT (delta_sol)
         x_sol2= xvent+SQRT (delta_sol)
         
         IF ( (xb(i-1) < x_sol1).AND. (x_sol1 < xb(i)) ) THEN
            n=n+1
            poly_p(n)%x=x_sol1
            poly_p(n)%y=yb(j)
            poly_p(n)%side=2
         ENDIF 
         
         IF ( (xb(i-1) < x_sol2) .AND. (x_sol2 < xb(i)) ) THEN
            n=n+1
            poly_p(n)%x=x_sol2
            poly_p(n)%y=yb(j)
            poly_p(n)%side=2
         ENDIF 
      ENDIF
      
      IF ( (xb(i)-xvent)**2 + (yb(j)-yvent)**2 <= vent_radius**2 ) THEN
         n=n+1
         poly_p(n)%x=xb(i)
         poly_p(n)%y=yb(j)
         poly_p(n)%side=0
      ENDIF
      
      delta_sol=vent_radius**2-(xb(i)-xvent)**2
      
      IF ( delta_sol > 0 ) THEN
         y_sol1= yvent+SQRT (delta_sol)
         y_sol2= yvent-SQRT (delta_sol)
         
         IF ( (yb(j-1) < y_sol1) .AND. (y_sol1 < yb(j)) ) THEN
            n=n+1
            poly_p(n)%x=xb(i)
            poly_p(n)%y=y_sol1
            poly_p(n)%side=3
         ENDIF 
         
         IF ( (yb(j-1) < y_sol2) .AND. ( y_sol2 < yb(j) ) ) THEN
            n=n+1
            poly_p(n)%x=xb(i)
            poly_p(n)%y=y_sol2
            poly_p(n)%side=3
         ENDIF 
      ENDIF
      
      IF ( (xb(i)-xvent)**2 + (yb(j-1)-yvent)**2 <= vent_radius**2 ) THEN
         n=n+1
         poly_p(n)%x=xb(i)
         poly_p(n)%y=yb(j-1)
         poly_p(n)%side=0
      ENDIF
      
      delta_sol=vent_radius**2-(yb(j-1)-yvent)**2
      
      IF ( delta_sol > 0 ) THEN
         x_sol1= xvent+SQRT (delta_sol)
         x_sol2= xvent-SQRT (delta_sol)
         IF ( (xb(i-1) < x_sol1).AND. (x_sol1 < xb(i)) ) THEN
            n=n+1
            poly_p(n)%x=x_sol1
            poly_p(n)%y=yb(j-1)
            poly_p(n)%side=4
         ENDIF 
         
         IF ( (xb(i-1) < x_sol2) .AND. (x_sol2 < xb(i)) ) THEN
            n=n+1
            poly_p(n)%x=x_sol2
            poly_p(n)%y=yb(j-1)
            poly_p(n)%side=4
         ENDIF 
      ENDIF
      
      DO l=1,n
         poly_p(n+l)%x=poly_p(l)%x
         poly_p(n+l)%y=poly_p(l)%y
         poly_p(n+l)%side=poly_p(l)%side
      ENDDO
      
      area_p=0.d0
      
      DO l=1,n-2
         area_p=area_p+areaT(poly_p(1)%x,poly_p(l+1)%x,poly_p(l+2)%x,   &
         poly_p(1)%y,poly_p(l+1)%y,poly_p(l+2)%y)
      ENDDO
      
      DO l=1,n
         IF  ( ( poly_p(l)%side*poly_p(l+1)%side > 0 ) .AND. &
             ( poly_p(l)%side /= poly_p(l+1)%side ) ) THEN
           
             area_p = area_p + areaSEC(poly_p(l)%x,poly_p(l+1)%x,  &
             poly_p(l)%y,poly_p(l+1)%y,xvent,yvent,vent_radius) &
             - areaT(poly_p(l)%x,poly_p(l+1)%x,xvent,poly_p(l)%y, &
             poly_p(l+1)%y,yvent)
             
          ENDIF
      ENDDO
      
      cell_fraction= area_p / abs( (xb(i-1)-xb(i))*(yb(j-1)-yb(j)) ) 
      
      RETURN
      END FUNCTION cell_fraction
!-----------------------------------------------------------------------
      REAL*8 FUNCTION areaT(x1,x2,x3,y1,y2,y3)
!
! ... evaluate the area of a triangle
!
      REAL*8 :: x1,x2,x3,y1,y2,y3
      REAL*8 :: p,a,b,c
      
      a=SQRT( (x1-x2)**2 + (y1-y2)**2 )
      b=SQRT( (x1-x3)**2 + (y1-y3)**2 )
      c=SQRT( (x3-x2)**2 + (y3-y2)**2 )
      p= (a+b+c) / 2.d0
      
      areaT= SQRT ( p * (p-a) * (p-b) * (p-c) )
      
      RETURN
      END FUNCTION areaT
!-----------------------------------------------------------------------
      REAL*8 FUNCTION areaSEC(x1,x2,y1,y2,xC,yC,R)
!
! ... evaluate the area of a sector of circle
!
      REAL*8 :: x1,x2,y1,y2,xC,yC,R
      REAL*8 :: alpha1,alpha2,alpha
      REAL*8 :: pi
      
      pi = 4.D0 * datan(1.D0)
      
      alpha1= datan2( (y1-yC)/R , (x1-xC)/R )
      alpha2= datan2( (y2-yC)/R , (x2-xC)/R )
      
      alpha= ABS (alpha1-alpha2)
      IF ( alpha > pi ) alpha = 2.D0*pi - alpha
      
      areaSEC = 0.5D0*alpha*R**2
      
      RETURN
      END FUNCTION areaSEC
!-----------------------------------------------------------------------
      END MODULE vent_conditions
!-----------------------------------------------------------------------
