!-----------------------------------------------------------------------
      MODULE dome_conditions
!
! ... Identify the dome cells
! ... Set the initial conditions within the dome cells
!
!-----------------------------------------------------------------------
      USE dimensions, ONLY: max_nsolid, max_ngas, nsolid, ngas
      USE io_files, ONLY: domeunit, domefile
      IMPLICIT NONE

      ! ... flags
      !
      INTEGER :: idome, idw, ishape, domegeom
      !
      ! ... 3 shells model parameters
      REAL*8 :: xdome, ydome, zdome
      REAL*8 :: volume_correction
      REAL*8 :: radius_1, radius_2, radius_3
      REAL*8 :: dome_radius
      REAL*8 :: mass_1, mass_2, mass_3
      REAL*8 :: rhos_1, rhos_2, rhos_3
      REAL*8 :: volume_1, volume_2, volume_3
      REAL*8 :: particle_fraction_1(max_nsolid)
      REAL*8 :: particle_fraction_2(max_nsolid)
      REAL*8 :: particle_fraction_3(max_nsolid)
      REAL*8 :: temperature_1, temperature_2, temperature_3
      REAL*8 :: overpressure_1, overpressure_2, overpressure_3
      !
      ! ... Parameters for Woods model
      REAL*8 :: gas_flux
      REAL*8 :: permeability
      REAL*8 :: dome_gasvisc
      REAL*8 :: conduit_radius
      !
      REAL*8 :: dome_ygc(max_ngas)

      TYPE icdome_cell
        INTEGER :: imesh
        INTEGER :: i
        INTEGER :: j
        INTEGER :: k
        REAL*8  :: radius
        REAL*8  :: angle
        REAL*8  :: pressure
        REAL*8  :: temperature
        REAL*8  :: sfraction(max_nsolid)
      END TYPE icdome_cell
      TYPE(icdome_cell), ALLOCATABLE :: dcell(:)

      INTEGER :: ndm, iid, jjd, kkd
      PUBLIC
      PRIVATE :: icdome_cell, dcell, ndm, iid, jjd, kkd
      
      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
! ... This routine locate the vent on a georeferenced mesh and set the
! ... cell flags onto the neighbours of the inlet
!
      SUBROUTINE locate_dome

      USE control_flags, ONLY: job_type, lpr
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE parallel, ONLY: mpime, root
      USE particles_constants, ONLY: rl

      IMPLICIT NONE
      
      INTEGER :: n
      REAL*8 :: pi
      REAL*8 :: dome_base
!
      pi = 4.D0 * ATAN(1.D0)
!
! ... Compute the dome actual radius, given the solid mass and concentrations,
! ... "volume correction" is the fraction of the sphere volume (=0.5 for a hemisphere)
! ... The three layers are representative of:
! ... 1. dome core 
! ... 2. transition layer or carapace
! ... 3. rock cap 
!
      volume_1 = 0.D0; volume_2 = 0.D0; volume_3 = 0.D0 
      !
      rhos_1 = SUM(particle_fraction_1(1:nsolid)*rl(1:nsolid))
      IF (rhos_1 /= 0.D0) volume_1 = mass_1 / rhos_1 / volume_correction
      !
      rhos_2 = SUM(particle_fraction_2(1:nsolid)*rl(1:nsolid))
      IF (rhos_2 /= 0.D0) volume_2 = mass_2 / rhos_2 / volume_correction
      !
      rhos_3 = SUM(particle_fraction_3(1:nsolid)*rl(1:nsolid))
      IF (rhos_3 /= 0.D0) volume_3 = mass_3 / rhos_3 / volume_correction
!
      IF (domegeom == 1) THEN
        ! ... For a spherical dome use the following radii
        radius_1 = ( 0.75D0 * volume_1 / pi )**(1.0/3.0)
        radius_2 = ( 0.75D0 * (volume_1 + volume_2) / pi )**(1.0/3.0)
        radius_3 = ( 0.75D0 * (volume_1 + volume_2 + volume_3) / pi )**(1.0/3.0)
      ELSE IF (domegeom == 2) THEN
        ! ... For a cylindrical dome use the following radii
        dome_base = pi * dome_radius**2
        radius_1 = volume_1 / dome_base
        radius_2 = (volume_1 + volume_2) / dome_base
        radius_3 = (volume_1 + volume_2 + volume_3) / dome_base
      END IF
!
      IF( job_type == JOB_TYPE_2D) THEN
              CALL locate_dome_2D
      ELSE IF( job_type == JOB_TYPE_3D) THEN
              CALL locate_dome_3D
      END IF
!
! ... Print out the dome coordinates
!
      IF( mpime == root ) THEN
        OPEN(UNIT=domeunit,FILE=domefile,STATUS='UNKNOWN')
        WRITE(domeunit,*) 'Report of dome initial conditions'
        WRITE(domeunit,*) 
        WRITE(domeunit,100) iid, jjd, kkd
        WRITE(domeunit,200) xdome, ydome, zdome
        WRITE(domeunit,400) radius_1, radius_2, radius_3
        WRITE(domeunit,*) 
        WRITE(domeunit,*) 'Number of dome cells: ', ndm
        WRITE(domeunit,*) 'Dome cells report'
        DO n = 1, ndm
          WRITE(domeunit,500) n, dcell(n)%imesh, dcell(n)%i, dcell(n)%j, dcell(n)%k, dcell(n)%radius, dcell(n)%angle
        END DO
100     FORMAT(1X,'dome center indices: ',3I5)
200     FORMAT(1X,'dome center coordinates: ',3(F12.2))
400     FORMAT(1X,'radii of the 3 shells: ',3(F12.2))
500     FORMAT(I6, I10, 3I5, F12.2, F12.6)
      END IF
!
      RETURN
      END SUBROUTINE locate_dome
!-----------------------------------------------------------------------
      SUBROUTINE locate_dome_2D

      USE atmospheric_conditions, ONLY: gravz
      USE dimensions, ONLY: nx, nz, ntot
      USE grid, ONLY: x, z, fl, xb, zb, dz, r, itc
      USE grid, ONLY: flag, fluid, dome_cell, slip_wall, noslip_wall, bl_cell
      USE volcano_topography, ONLY: itp, ord
      USE particles_constants, ONLY: rl

      IMPLICIT NONE
      
      INTEGER :: i, j, k, ijk, n, is
      REAL*8 :: distance, distance2, rp
!
! ... In cylindrical coordinates, the dome center is located on the axis.
! ... Otherwise its coordinates '(xdome)' are taken from the input.
!
      jjd = 1
      IF (itc >= 1 .OR. iid==1) THEN
        iid = 2
      ELSE IF (itc == 0) THEN
        DO i = 1, nx
          IF (x(i) <= xdome) iid = i
        END DO
      END IF
      xdome = x(iid)
!
! ... Define the quota of the volcanic dome
! ... (considering the topography). If the topography
! ... has been red from a file, 'kkd' has already been set
! ... in the 'volcano_topography' module.
!
      IF( itp >= 1 ) THEN
        kkd = ord(iid)
      ELSE
        DO k= 1, nz
          ijk = iid + (k-1) * nx 
          IF (fl(ijk) == slip_wall .OR. fl(ijk) == noslip_wall) kkd = k
        END DO
      END IF
      zdome = z(kkd)
!
! ... Count the cells of the dome and set the 'dome-cell' flags
!
      ndm = 0
      DO k = 1, nz
        DO i = 1, nx
          ijk = i + (k-1) * nx
          !
          IF (domegeom == 1) THEN
            ! ... Spherical dome
            distance2 = (x(i)-x(iid))**2 + (z(k)-z(kkd))**2
            IF (distance2 <= radius_3**2 .AND. (fl(ijk)==fluid .OR. fl(ijk)==bl_cell)) THEN
              ndm = ndm + 1
              fl(ijk) = dome_cell
            END IF
          ELSE IF (domegeom == 2) THEN
            ! ... Cylindrical dome
            distance2 = DABS(x(i)-x(iid))
            IF ((distance2 < dome_radius .AND. z(k)-z(kkd) < radius_3) .AND. &
                (fl(ijk)==fluid .OR. fl(ijk)==bl_cell)) THEN
              ndm = ndm + 1
              fl(ijk) = dome_cell
            END IF
          END IF
        END DO
      END DO
!      
! ... Allocate the data-type representing dome cells
! ... Allocate the array for the dome pressure profile
!
      ALLOCATE(dcell(ndm))
      dcell(:)%imesh = 0
      dcell(:)%i = 1
      dcell(:)%j = 1
      dcell(:)%k = 1
      dcell(:)%radius = 0.D0
      dcell(:)%angle = 0.D0
      dcell(:)%pressure = 0.D0
      dcell(:)%temperature = 0.D0
      DO is = 1, max_nsolid
        dcell(:)%sfraction(is) = 0.D0
      END DO
!
! ... Map the dome cells 
!
      n = 0 
      DO k = 1, nz
        DO i = 1, nx
          ijk = i + (k-1) * nx
          IF (fl(ijk) == dome_cell) THEN
                  distance = DSQRT( (x(i)-x(iid))**2 + (z(k)-z(kkd))**2 )
                  n = n + 1
                  dcell(n)%imesh = ijk
                  dcell(n)%i = i
                  dcell(n)%k = k
                  dcell(n)%radius = distance
                  !
                  IF (distance > 0.D0 ) THEN
                          dcell(n)%angle  = ATAN2((z(k)-z(kkd)),(x(i)-x(iid)))
                  ELSE
                          dcell(n)%angle = 0.D0
                  END IF
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE locate_dome_2D
!----------------------------------------------------------------------
! ... This routine locate the vent on a georeferenced mesh and set the
! ... cell flags onto the neighbours of the inlet
!
      SUBROUTINE locate_dome_3D

      USE atmospheric_conditions, ONLY: gravz
      USE dimensions, ONLY: nx, ny, nz, ntot
      USE grid, ONLY: x, y, z, fl, xb, yb, zb, dz, r, itc
      USE grid, ONLY: flag, fluid, dome_cell, slip_wall, noslip_wall, bl_cell
      USE volcano_topography, ONLY: itp, ord2d
      USE volcano_topography, ONLY: flatten_dem_dome
      USE particles_constants, ONLY: rl

      IMPLICIT NONE
      
      INTEGER :: i, j, k, ijk, n, is
      REAL*8 :: distance, distance2, rp
!
! ... The dome coordinates '(xdome, ydome)' are always taken from the input.
!
      DO i = 1, nx
        IF (x(i) <= xdome) iid = i
      END DO
      DO j = 1, ny
        IF (y(j) <= ydome) jjd = j
      END DO
      xdome = x(iid)
      ydome = y(jjd)
!
! ... Define the quota of the volcanic dome
! ... (considering the topography). If the topography
! ... has been red from a file, 'kkd' has already been set
! ... in the 'volcano_topography' module.
!
      IF( itp >= 1 ) THEN
        kkd = ord2d(iid,jjd)
      ELSE
        DO k= 1, nz
          ijk = iid + (jjd-1) * nx + (k-1) * nx * ny
          IF (fl(ijk) == slip_wall .OR. fl(ijk) == noslip_wall) kkd = k
        END DO
      END IF
!
! ... If the dome center is lower than the topographic elevation
! ... the topography is "excavated". Otherwise the center is placed
! ... on the topography
!
      IF (zdome < zb(kkd)) THEN
              DO k = 1, nz
                IF (zb(k) <= zdome) kkd = k
              END DO
              !
              IF (itp >=1 ) THEN
                IF (domegeom == 1) THEN
                  ! ... Spherical dome
                  CALL flatten_dem_dome(xdome,ydome,radius_3,kkd)
                ELSE IF (domegeom == 2) THEN
                  ! ... Cylindrical dome
                  CALL flatten_dem_dome(xdome,ydome,dome_radius,kkd)
                END IF
              END IF
      ELSE
              zdome = z(kkd)
      END IF
!
! ... Count the cells of the dome
!
      ndm = 0
      DO k = 1, nz
        DO j = 1, ny
          DO i = 1, nx
            ijk = i + (j-1) * nx + (k-1) * nx * ny
            !
            IF (domegeom == 1) THEN
              ! ... Spherical dome
              distance2 = (x(i)-x(iid))**2 + (y(j)-y(jjd))**2 + (z(k)-z(kkd))**2
              IF (distance2 <= radius_3**2 .AND. (fl(ijk) == fluid .OR. fl(ijk)==bl_cell)) THEN
                    ndm = ndm + 1
                    fl(ijk) = dome_cell
              END IF
            ELSEIF (domegeom == 2) THEN
              ! ... Cylindrical dome
              distance2 = (x(i)-x(iid))**2 + (y(j)-y(jjd))**2
              IF ( (distance2 <= dome_radius**2 .AND. (z(k)-z(kkd)) <= radius_3) .AND.&
                  (fl(ijk)==fluid .OR. fl(ijk)==bl_cell)) THEN
                    ndm = ndm + 1
                    fl(ijk) = dome_cell
              END IF
            END IF
          END DO
        END DO
      END DO
!      
! ... Allocate the data-type representing dome cells
!
      ALLOCATE(dcell(ndm))
      dcell(:)%imesh = 0
      dcell(:)%i = 1
      dcell(:)%j = 1
      dcell(:)%k = 1
      dcell(:)%radius = 0.D0
      dcell(:)%angle = 0.D0
      dcell(:)%pressure = 0.D0
      dcell(:)%temperature = 0.D0
      DO is = 1, max_nsolid
        dcell(:)%sfraction(is) = 0.D0
      END DO
!
! ... Map the dome cells 
!
      n = 0
      DO k = 1, nz
        DO j = 1, ny
          DO i = 1, nx
            ijk = i + (j-1) * nx + (k-1) * nx * ny
            IF (fl(ijk) == dome_cell) THEN
                    distance = DSQRT( (x(i)-x(iid))**2 + (y(j)-y(jjd))**2 + (z(k)-z(kkd))**2 )
                    rp = DSQRT( (x(i)-x(iid))**2 + (y(j)-y(jjd))**2 )
                    n = n + 1
                    dcell(n)%imesh = ijk
                    dcell(n)%i = i
                    dcell(n)%j = j
                    dcell(n)%k = k
                    dcell(n)%radius = distance
                    IF (distance > 0.D0) THEN
                          dcell(n)%angle  = ATAN2((z(k)-z(kkd)), rp)
                    ELSE
                          dcell(n)%angle = 0.D0
                    END IF
            END IF
          END DO
        END DO
      END DO
!      
      RETURN
      END SUBROUTINE locate_dome_3D
!-----------------------------------------------------------------------
      SUBROUTINE dome_reshape
      USE domain_mapping, ONLY: ncint, meshinds, data_exchange
      USE domain_mapping, ONLY: myijk, myinds, set_myinds
      USE grid, ONLY: flag, dome_cell, noslip_wall
      USE io_files, ONLY: testunit
      USE time_parameters, ONLY: itd
      IMPLICIT NONE
      INTEGER :: imesh, ijk, i, j, k, n, counter
      LOGICAL :: condition
! 
      condition = .TRUE.
      counter = 0
      mesh_loop: DO ijk = 1, ncint      
          CALL meshinds(ijk,imesh,i,j,k)
          !
          ! ... Loop over the dome cells to find the
          ! ... dome-cell index 
          !
          search_n: DO n = 1, ndm
            IF (dcell(n)%imesh == imesh) THEN
              ! ... Condition of MSH wedged dome
              condition = ( 3*(dcell(n)%j-jjd)<(dcell(n)%k-kkd) )
              !condition = ( (dcell(n)%k-kkd) > 20 ) ! this works fine !
              !condition = ( dcell(n)%j-jjd < 0 .AND. dcell(n)%radius >= radius_2 )
              IF (condition .OR. itd>1) THEN
                flag(ijk) = noslip_wall
                counter = counter + 1
              END IF
              EXIT search_n
            END IF
          END DO search_n
      END DO mesh_loop
      WRITE(testunit,*) counter, 'dome cells blocked'
!
      CALL data_exchange(flag)
      CALL set_myinds(myinds, myijk)
!
      RETURN
      END SUBROUTINE dome_reshape
!-----------------------------------------------------------------------
      SUBROUTINE set_domec
!
! ... Compute the initial conditions for a semi-spherical dome
!
      USE array_filters, ONLY: interp
      USE atmospheric_conditions, ONLY: p_atm, gravz
      USE control_flags, ONLY: job_type, lpr
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nsolid, ngas, nx, ny, nz
      USE domain_mapping, ONLY: ncint, meshinds, data_exchange
      USE environment, ONLY: cpclock
      USE eos_gas, ONLY: ygc
      USE gas_constants, ONLY: gas_type, rgas
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_temperature, ONLY: tg, ts
      USE gas_solid_velocity, ONLY: ug, wg, vg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE grid, ONLY: flag, x, y, z, dome_cell, immb_cell, filled_cell_2
      USE grid, ONLY: dz
      USE io_files, ONLY: testunit
      USE parallel, ONLY: mpime, root
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: ep, p
      USE set_indexes, ONLY: ijkp, first_subscr
      IMPLICIT NONE

      INTEGER, ALLOCATABLE, DIMENSION(:) :: domeindex
      REAL*8, ALLOCATABLE :: sdensity(:), dpressure(:)
      REAL*8 :: ygcsum, ra, pi, psi, raddo, erre
      REAL*8 :: beta, p_hydro, p_ext
      REAL*8 :: grad, fact
      INTEGER :: ijk, imesh, i,j,k, is, ig, n, counter
      INTEGER :: kk, nn, ka
!      
! ... Set constant for the computation of dome pressure
!
      pi = 4.D0 * ATAN(1.D0)
      psi = 2.D0 * pi
!
! ... Allocate the arrays for dome cells
!
      ALLOCATE(sdensity(ndm),dpressure(ndm))
      sdensity = 0.D0 
      dpressure = 0.D0
!
! ... Allocate the array of dome indices
!
      ALLOCATE(domeindex(ncint))
      domeindex = 0
!
! ... Useful constants for the permeable dome model
! ... [Woods et al., 2002]
!
      IF (idome == 1) THEN
        p_ext = p_atm(kkd)
        erre = rgas / 18.D0 * 1.D3
        beta = 2.D0 * dome_gasvisc * gas_flux * erre * temperature_1
        beta = beta / ( p_ext**2 * permeability * psi * radius_2 )
      END IF
!
! ... Set initial conditions in the dome cells
! ... (Immersed Boundaries are NOT considered within the dome!)
!
      mesh_loop_1: DO ijk = 1, ncint      
        CALL first_subscr(ijk)
        IF(flag(ijk) == dome_cell .OR. &
          (flag(ijk) == immb_cell .AND. flag(ijkp)==dome_cell) .OR. &
          (flag(ijk) == filled_cell_2 .AND. flag(ijkp)==dome_cell)) THEN
          CALL meshinds(ijk,imesh,i,j,k)
          !
          ! ... Loop over the dome cells to find the
          ! ... dome-cell index 
          !
          n = 0
          search_n: DO counter = 1, ndm
            IF (dcell(counter)%imesh == imesh) THEN
              domeindex(ijk) = counter
              n = counter
              EXIT search_n
            END IF
          END DO search_n

          IF (n/=0) THEN
            ! ... Set the temperature, gas pressure and particle fractions
            ! ... in every cell of the dome
            !
            IF (domegeom == 1) THEN
              ! ... Spherical dome
              ra = dcell(n)%radius
            ELSE IF (domegeom == 2) THEN
              ! ... Cylindrical dome
              ra = z(k) - zdome
            END IF
            
            ! ... Apply different properties for each shells 
            !
            IF (ra <= radius_1) THEN
                    dcell(n)%temperature = temperature_1
                    dcell(n)%pressure = overpressure_1
                    dcell(n)%sfraction(1:nsolid) = particle_fraction_1(1:nsolid)
            ELSE IF (ra <= radius_2 .AND. ra > radius_1) THEN
                    !
                    ! constant values
                    !
                    ! dcell(n)%temperature = temperature_2
                    ! dcell(n)%pressure = overpressure_2
                    ! dcell(n)%sfraction(1:nsolid) = particle_fraction_2(1:nsolid)
                    !
                    ! linear interpolation
                    !
                    fact = (ra-radius_1)/(radius_2-radius_1)
                    dcell(n)%temperature = temperature_1 + (temperature_3-temperature_1)*fact
                    dcell(n)%pressure = overpressure_1 + (overpressure_3-overpressure_1)*fact
                    dcell(n)%sfraction(1:nsolid) = particle_fraction_1(1:nsolid) + &
                      (particle_fraction_3(1:nsolid)-particle_fraction_1(1:nsolid))*fact
            ELSE IF (ra <= radius_3 .AND. ra > radius_2) THEN
                    dcell(n)%temperature = temperature_3
                    dcell(n)%pressure = overpressure_3
                    dcell(n)%sfraction(1:nsolid) = particle_fraction_3(1:nsolid)
            END IF
            sdensity(n) = SUM(dcell(n)%sfraction(1:nsolid)*rl(1:nsolid))
            !
          END IF
        END IF
      END DO mesh_loop_1
!
! ... All processors must know the density distribution to
! ... compute the hydrostatic pressure
!
      CALL parallel_sum_real(sdensity,ndm)
!
      mesh_loop_2: DO ijk = 1, ncint      
        !
        ! ... Set up of flow variables in the dome cells
        ! ... and compute the hydrostatic pressure (for solids) if prescribed
        !
        CALL first_subscr(ijk)
        IF(flag(ijk) == dome_cell .OR. &
          (flag(ijk) == immb_cell .AND. flag(ijkp)==dome_cell) .OR. &
          (flag(ijk) == filled_cell_2 .AND. flag(ijkp)==dome_cell)) THEN
          CALL meshinds(ijk,imesh,i,j,k)
          ra = DSQRT((x(i)-x(iid))**2+(y(j)-y(jjd))**2+(z(k)-z(kkd))**2)
          !
          ! ... Loop over the dome cells to find the
          ! ... dome-cell index 
          !
          n = domeindex(ijk)
            ! ... Set the initial conditions, as
            ! ... specified in the input file on 
            ! ... all cells enclosing the vent
            !
          IF (n/=0) THEN
            !
            ug(ijk) = 0.D0
            IF (job_type == JOB_TYPE_3D) vg(ijk) = 0.D0 
            wg(ijk) = 0.D0
            !
            tg(ijk) = dcell(n)%temperature
            !
            ep(ijk) = 1.D0
            DO is = 1, nsolid
              ep(ijk) = ep(ijk) - dcell(n)%sfraction(is)
              rlk(ijk,is) = dcell(n)%sfraction(is)*rl(is)
              !
              us(ijk,is)  = 0.D0
              IF (job_type == JOB_TYPE_3D) vs(ijk,is)  = 0.D0
              ws(ijk,is) = 0.D0
              !
              ! homogenous temperature applied in each shell/layer
              ts(ijk,is)  = dcell(n)%temperature
            END DO
            !
            ! ... Set the dome pressure
            IF (idome == 1) THEN
              ! ... Woods' (2002) profile
              p(ijk)  = p_dome(ra,p_ext,beta)
            ELSE IF (idome == 2) THEN
              ! ... Constant overpressure
              p(ijk) = p(ijk) + dcell(n)%pressure
            END IF
            !
            ! ... Add the hydrostatic pressure due to dome mass
            !  
            p_hydro = 0.D0
            IF (idw >= 1) THEN
              search_dome: DO kk = k, nz
                IF (job_type == JOB_TYPE_2D) THEN
                  imesh = i + (kk-1) * nx
                ELSE IF (job_type == JOB_TYPE_3D) THEN
                  imesh = i + (j-1) * nx + (kk-1) * nx * ny
                END IF
                !
                ! ... Loop again over the dome cells to
                ! ... find the index
                nn = 0
                DO counter = 1, ndm
                  IF (dcell(counter)%imesh == imesh) nn = counter
                END DO 
                !
                IF (nn /= 0) THEN
                  p_hydro = p_hydro - sdensity(nn)*gravz*dz(kk)
                ELSE
                    EXIT search_dome
                END IF
              END DO search_dome
              !
              p(ijk)  = p(ijk) + p_hydro
            END IF
            !
            DO ig = 1, ngas
              ygc(ijk,ig) = dome_ygc(gas_type(ig))
            END DO
            ! ... check gas components closure relation
            !
            ygcsum = SUM(ygc(ijk,:))
            IF ( ygcsum /= 1.D0 ) THEN
              ygc(ijk,ngas) = 1.D0 - SUM( ygc(ijk,1:ngas-1) )
            END IF
!
! ... update the local dome cell pressure 
!
            dcell(n)%pressure = p(ijk)
            dpressure(n) = dcell(n)%pressure
          END IF
!          
        END IF
      END DO mesh_loop_2
!
! ... Update the dome cell pressure in all processors
!
      CALL parallel_sum_real(dpressure, ndm)
      DO n = 1, ndm
        dcell(n)%pressure = dpressure(n)
      END DO
!
! ... Report dome initial conditions
!
      IF (mpime == root) THEN
              WRITE(domeunit,*) 
              IF (idw == 1) &
                WRITE(domeunit,*) 'Hydrostatic pressure is computed'
              !
              IF (idome == 1) THEN
                WRITE(domeunit,*) 
                WRITE(domeunit,*) 'Woods radial pressure profile'
                raddo = 0.D0
                !
                ! rocks mixed into dome particles 
                DO WHILE (raddo <= radius_3)
                  WRITE(domeunit,*) raddo, p_dome(raddo,p_ext,beta)
                  raddo = raddo + 1.D0
                END DO
              ELSE IF (idome == 2) THEN
                WRITE(domeunit,*) 'Constant overpressure (shell 1) = ', overpressure_1
                WRITE(domeunit,*) 'Constant overpressure (shell 2) = ', overpressure_2
                WRITE(domeunit,*) 'Constant overpressure (shell 3) = ', overpressure_3
              END IF
              !
              WRITE(domeunit,*) 
              WRITE(domeunit,*) 'Dome radial pressure profile'
              DO n = 1, ndm
                IF (dcell(n)%i == iid .AND. dcell(n)%j == jjd) THEN
                  WRITE(domeunit,*) dcell(n)%radius, dcell(n)%pressure
                END IF
              END DO
      END IF
 100  FORMAT(F12.2, G12.4)
      CALL compute_dome_mass_energy
!
! ... If .TRUE., modify dome shape 
! 
      IF (.FALSE.) THEN
        CALL dome_reshape
        CALL compute_dome_mass_energy
      END IF
!
      DEALLOCATE(sdensity)
      DEALLOCATE(dpressure)
      DEALLOCATE(domeindex)
!
      RETURN
      END SUBROUTINE set_domec
!-----------------------------------------------------------------------
      SUBROUTINE compute_dome_mass_energy
      USE atmospheric_conditions, ONLY: p_atm
      USE control_flags, ONLY: job_type, lpr
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nx, ny, nz
      USE domain_mapping, ONLY: ncint, meshinds
      USE gas_constants, ONLY: rgas
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_temperature, ONLY: tg, ts
      USE grid, ONLY: x, y, z, xb, yb, zb, r
      USE grid, ONLY: flag, dome_cell, immb_cell
      USE parallel, ONLY: mpime, root
      USE particles_constants, ONLY: rl
      USE pressure_epsilon, ONLY: p
      USE set_indexes, ONLY: first_subscr, ijkp
      IMPLICIT NONE
      INTEGER :: ijk, i, j, k, imesh, counter, n
      REAL*8 :: cgam, gammawater, muwater
      REAL*8 :: pi, twopi
      REAL*8 :: isoe, adie, mgd, mpd, vold, dvol
      REAL*8 :: mg1, mg2, mg3, mp1, mp2, mp3, mg, mp
      REAL*8 :: soldfr, voidfr, fact, distance

      gammawater = 1.33D0
      muwater    = 18.D-3
      cgam = ( gammawater - 1.D0 ) / gammawater
      pi = 4.D0 * ATAN(1.D0)
      twopi = 2.D0 * pi
!
      vold = 0.D0
      mgd = 0.D0
      mpd = 0.D0
      adie = 0.D0
      isoe = 0.D0
      mg1 = 0.D0; mg2 = 0.D0; mg3 = 0.D0 
      mp1 = 0.D0; mp2 = 0.D0; mp3 = 0.D0 
!
      mesh_loop: DO ijk = 1, ncint      
        CALL first_subscr(ijk)
        IF(flag(ijk) == dome_cell .OR. &
          (flag(ijk) == immb_cell .AND. flag(ijkp) == dome_cell)) THEN
          CALL meshinds(ijk,imesh,i,j,k)
          !
          ! ... Loop over the dome cells to find the
          ! ... cell index (must be optimized!!!)
          !
          n = 0
          DO counter = 1, ndm
            IF (dcell(counter)%imesh == imesh) n = counter
          END DO 
          IF (n/=0) THEN
            IF (domegeom == 1) THEN
              distance = dcell(n)%radius
            ELSE IF (domegeom == 2) THEN
              distance = z(k) - z(kkd)
            END IF
            fact = rgas * dcell(n)%temperature / muwater 
            soldfr = SUM( dcell(n)%sfraction(1:nsolid) )
            voidfr = 1.D0 - soldfr
            !
            ! ... volume of a cell
            IF (job_type == JOB_TYPE_2D) THEN
              dvol = twopi * r(i) * (xb(i)-xb(i-1))*(zb(k)-zb(k-1))
            ELSE IF (job_type == JOB_TYPE_3D) THEN
              dvol = (xb(i)-xb(i-1))*(yb(j)-yb(j-1))*(zb(k)-zb(k-1))
            END IF
            !
            ! ... specific energy 
            ! ... [ Neri et al. 1998, JVGR ]
            ! ... [ Wilson et al. 1980, Geophys. J. R. Astron. Soc. ]
            isoe = isoe + voidfr * p(ijk) * log(p(ijk)/p_atm(k)) * dvol
            adie = adie + voidfr * p(ijk) * (1.D0 - (p_atm(k)/p(ijk))**cgam) * dvol
            !
            ! ... mass per unit volume (cell)
            mg = p(ijk) * dvol / fact * voidfr
            mp = SUM(rlk(ijk,1:nsolid)) * dvol
            !
            IF (distance <= radius_1) THEN
              mg1 = mg1 + mg
              mp1 = mp1 + mp
            ELSE IF (distance <= radius_2 .AND. distance > radius_1) THEN
              mg2 = mg2 + mg
              mp2 = mp2 + mp
            ELSE IF (distance <= radius_3 .AND. distance > radius_2) THEN
              mg3 = mg3 + mg
              mp3 = mp3 + mp
            END IF
            !
            ! ... volume
            vold = vold + dvol
          ELSE 
            distance = 0.D0
            fact = rgas * tg(ijk) / muwater 
            soldfr = SUM( rlk(ijk,1:nsolid)/rl(1:nsolid) )
            voidfr = 1.D0 - soldfr
            !
            ! ... volume of a cell
            IF (job_type == JOB_TYPE_2D) THEN
              dvol = twopi * r(i) * (xb(i)-xb(i-1))*(zb(k)-zb(k-1))
            ELSE IF (job_type == JOB_TYPE_3D) THEN
              dvol = (xb(i)-xb(i-1))*(yb(j)-yb(j-1))*(zb(k)-zb(k-1))
            END IF
            !
            ! ... specific energy 
            ! ... [ Neri et al. 1998, JVGR ]
            ! ... [ Wilson et al. 1980, Geophys. J. R. Astron. Soc. ]
            isoe = isoe + voidfr * p(ijk) * log(p(ijk)/p_atm(k)) * dvol
            adie = adie + voidfr * p(ijk) * (1.D0 - (p_atm(k)/p(ijk))**cgam) * dvol
            !
            ! ... mass per unit volume (cell)
            mg = p(ijk) * dvol / fact * voidfr
            mp = SUM(rlk(ijk,1:nsolid)) * dvol
            WRITE(*,*) 'proc ', mpime, 'cell: ', ijk, fact, soldfr, voidfr, mg, mp
            !
              mg1 = 0.D0
              mp1 = 0.D0
              mg2 = 0.D0
              mp2 = 0.D0
              mg3 = mg3 + mg
              mp3 = mp3 + mp
            !
            ! ... volume
            vold = vold + dvol
          END IF
        END IF
      END DO mesh_loop
!
      CALL parallel_sum_real(isoe,1)
      CALL parallel_sum_real(adie,1)
      CALL parallel_sum_real(mg1,1)
      CALL parallel_sum_real(mp1,1)
      CALL parallel_sum_real(mg2,1)
      CALL parallel_sum_real(mp2,1)
      CALL parallel_sum_real(mg3,1)
      CALL parallel_sum_real(mp3,1)
      CALL parallel_sum_real(vold,1)
      !
      ! ... solid and gas mass
      !
      mgd = mg1 + mg2 + mg3
      mpd = mp1 + mp2 + mp3
      ! ... 
      isoe = isoe / (mpd + mgd)
      adie = adie / cgam / (mpd + mgd)
!
! ... Print out the dome coordinates
!
      IF(mpime == root ) THEN
              WRITE(domeunit,*)
              WRITE(domeunit,*) 'Mass and Energy stored in the dome'
              WRITE(domeunit,100) vold 
              WRITE(domeunit,198) mg1, mg2, mg3
              WRITE(domeunit,199) mp1, mp2, mp3
              WRITE(domeunit,200) mgd
              WRITE(domeunit,201) mpd
              WRITE(domeunit,300) isoe
              WRITE(domeunit,400) adie
              WRITE(domeunit,*)
100     FORMAT(1X,'Dome Volume: ',F18.4)
198     FORMAT(1X,'Gas   mass [Kg]: ',3F18.4)
199     FORMAT(1X,'Solid   mass [Kg]: ',3F18.4)
200     FORMAT(1X,'Total gas   mass [Kg]: ',F18.4)
201     FORMAT(1X,'Total solid mass [Kg]: ',F18.4)
300     FORMAT(1X,'Specific Isothermal expansion energy [J/Kg]: ',F18.4)
400     FORMAT(1X,'Specific Adiabatic  expansion energy [J/Kg]: ',F18.4)
      END IF

      RETURN
      END SUBROUTINE compute_dome_mass_energy
!-----------------------------------------------------------------------
      REAL*8 FUNCTION dh(r, a)
      !
      ! ... Distance of a point in the dome from the top surface
      !
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: r, a
      REAL*8 :: a1

      a1 = ACOS(r/radius_3 * COS(a))
      dh = radius_3 * SIN(a1) - r * SIN(a)

      RETURN
      END FUNCTION dh
!-----------------------------------------------------------------------
      REAL*8 FUNCTION p_dome(er,pa,beta)
      !
      ! ... Dome pressurization model (Woods et al., 2002)
      !
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: er, pa, beta
      REAL*8 :: cr, fact

      cr = MAX(er, conduit_radius)
      IF (cr <= radius_2) THEN
              fact = radius_2 / cr - 1.D0
      ELSE
              fact = 0.D0
      END IF
      p_dome = pa * DSQRT(1.D0 + beta * fact)

      RETURN
      END FUNCTION p_dome
!-----------------------------------------------------------------------
      END MODULE dome_conditions
!-----------------------------------------------------------------------
