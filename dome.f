!-----------------------------------------------------------------------
      MODULE dome_conditions
!-----------------------------------------------------------------------
      USE dimensions, ONLY: max_nsolid, max_ngas, nsolid, ngas
      USE io_files, ONLY: logunit
      IMPLICIT NONE

      ! ... flags
      !
      INTEGER :: idome, idw

      REAL*8 :: xdome, ydome, zdome, dome_volume, dome_radius
      REAL*8 :: temperature, overpressure
      REAL*8 :: particle_fraction(max_nsolid)
      !
      REAL*8 :: dome_ygc(max_ngas)

      TYPE icdome_cell
        INTEGER :: imesh
        REAL*8  :: radius
        REAL*8  :: angle
        REAL*8  :: p_hydro
      END TYPE icdome_cell

      TYPE(icdome_cell), ALLOCATABLE :: dcell(:)

      INTEGER :: ndm
      REAL*8 :: gas_flux
      REAL*8 :: permeability
      REAL*8 :: dome_gasvisc
      REAL*8 :: conduit_radius
      PUBLIC
      
      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
! ... This routine locate the vent on a georeferenced mesh and set the
! ... cell flags onto the neighbours of the inlet
!
      SUBROUTINE locate_dome

      USE atmospheric_conditions, ONLY: gravz
      USE control_flags, ONLY: job_type, lpr
      USE dimensions, ONLY: nx, ny, nz, ntot
      USE grid, ONLY: x, y, z, fl, xb, yb, zb, dz, r, itc
      USE grid, ONLY: iv, jv, kv, grigen
      USE grid, ONLY: flag, fluid, dome_cell, slip_wall, noslip_wall
      USE volcano_topography, ONLY: itp, ord2d
      USE volcano_topography, ONLY: flatten_dem
      USE parallel, ONLY: mpime, root
      USE particles_constants, ONLY: rl

      IMPLICIT NONE
      
      INTEGER :: i, j, k, ijk, n
      REAL*8 :: distance, distance2, pi, twopi, rp
      REAL*8 :: rlks
!
! ... Compute the dome radius, given the volume.
! ... The dome shape is half-a-sphere. 
!
      pi = 4.D0 * ATAN(1.D0)
      twopi = 2.D0 * pi
      dome_radius = ( 1.5D0 * dome_volume / pi )**(1.0/3.0)
      rlks = SUM( particle_fraction(1:nsolid)*rl(1:nsolid) ) 
!
      IF( job_type == '2D') THEN
!
! ... If the mesh has been generated by the code,
! ... the dome center is located where the mesh is more refined.
! ... Otherwise its coordinates '(xdome)' are taken from the input.
!
        IF (itc >= 1 .OR. iv==1) iv = 2
        IF (grigen >= 1) THEN
          xdome = x(iv)
        ELSE
          DO i = 1, nx
            IF (x(i) <= xdome) iv = i
          END DO
          xdome = x(iv)
        END IF
!
! ... Define the quota of the volcanic dome
! ... (considering the topography). If the topography
! ... has been red from a file, 'kv' has already been set
! ... in the 'volcano_topography' module.
!
        IF( itp == 0 ) THEN
          DO k= 1, nz
            ijk = iv + (k-1) * nx 
            IF (fl(ijk) == slip_wall .OR. fl(ijk) == noslip_wall) kv = k
          END DO
        END IF
        zdome = z(kv)
!
! ... Count the cells of the dome and set the 'dome-cell' flags
!
        ndm = 0
        DO k = 1, nz
          DO i = 1, nx
            ijk = i + (k-1) * nx
            distance2 = (x(i)-x(iv))**2 + (z(k)-z(kv))**2
            IF (distance2 <= dome_radius**2 .AND. fl(ijk) == fluid ) THEN
                    ndm = ndm + 1
                    fl(ijk) = dome_cell
            END IF
          END DO
        END DO
!      
! ... Allocate the data-type representing dome cells
!
        ALLOCATE(dcell(ndm))
!
! ... Map the dome cells
!
        n = 0
        DO k = 1, nz
          DO i = 1, nx
            ijk = i + (k-1) * nx
            IF (fl(ijk) == dome_cell) THEN
                    distance = DSQRT( (x(i)-x(iv))**2 + (z(k)-z(kv))**2 )
                    n = n + 1
                    dcell(n)%imesh = ijk
                    dcell(n)%radius = distance
                    IF (distance > 0.D0 ) THEN
                            dcell(n)%angle  = ATAN2((z(k)-z(kv)),(x(i)-x(iv)))
                    ELSE
                            dcell(n)%angle = 0.D0
                    END IF
                    dcell(n)%p_hydro  = rlks * DABS(gravz) * dh(distance,dcell(n)%angle)
            END IF
          END DO
        END DO
!      
      ELSE IF (job_type == '3D') THEN
!
! ... If the mesh has been generated by the code,
! ... the dome center is located where the mesh is more refined.
! ... Otherwise its coordinates '(xdome, ydome)' are taken from the input.
!
        IF (grigen >= 1) THEN
          xdome = x(iv)
          ydome = y(jv)
        ELSE 
          DO i = 1, nx
            IF (x(i) <= xdome) iv = i
          END DO
          DO j = 1, ny
            IF (y(j) <= ydome) jv = j
          END DO
          xdome = x(iv)
          ydome = y(jv)
        END IF
!
! ... Define the quota of the volcanic dome
! ... (considering the topography). If the topography
! ... has been red from a file, 'kv' has already been set
! ... in the 'volcano_topography' module.
!
        IF( itp < 1 ) THEN
          DO k= 1, nz
            ijk = iv + (jv-1) * nx + (k-1) * nx * ny
            IF (fl(ijk) == slip_wall .OR. fl(ijk) == noslip_wall) kv = k
          END DO
        END IF
        zdome = z(kv)
!
! ... Count the cells of the dome
!
        ndm = 0
        DO k = 1, nz
          DO j = 1, ny
            DO i = 1, nx
              ijk = i + (j-1) * nx + (k-1) * nx * ny
              distance2 = (x(i)-x(iv))**2 + (y(j)-y(jv))**2 + (z(k)-z(kv))**2
              IF (distance2 <= dome_radius**2 .AND. fl(ijk) == fluid ) THEN
                      ndm = ndm + 1
                      fl(ijk) = dome_cell
              END IF
            END DO
          END DO
        END DO
!      
! ... Allocate the data-type representing dome cells
!
        ALLOCATE(dcell(ndm))
!
! ... Map the dome cells
!
        n = 0
        DO k = 1, nz
          DO j = 1, ny
            DO i = 1, nx
              ijk = i + (j-1) * nx + (k-1) * nx * ny
              IF (fl(ijk) == dome_cell) THEN
                      distance = DSQRT( (x(i)-x(iv))**2 + (y(j)-y(jv))**2 + (z(k)-z(kv))**2 )
                      rp = DSQRT( (x(i)-x(iv))**2 + (y(j)-y(jv))**2 )
                      n = n + 1
                      dcell(n)%imesh = ijk
                      dcell(n)%radius = distance
                      IF (distance > 0.D0) THEN
                            dcell(n)%angle  = ATAN2((z(k)-z(kv)), rp)
                      ELSE
                            dcell(n)%angle = 0.D0
                      END IF
                      dcell(n)%p_hydro  = rlks * DABS(gravz) * dh(distance, dcell(n)%angle)
              END IF
            END DO
          END DO
        END DO
!      
      END IF
!
! ... Print out the dome coordinates
!
      IF( lpr > 0 .AND. mpime == root ) THEN
        WRITE(logunit,100) iv, jv, kv
        WRITE(logunit,200) xdome, ydome, zdome
        WRITE(logunit,400) dome_radius
100     FORMAT(1X,'dome center: ',3I5)
200     FORMAT(1X,'dome center coordinates: ',3(F12.2))
300     FORMAT(1X,'dome volume: ',F12.2)
400     FORMAT(1X,'dome radius: ',F12.2)
      END IF
!
      RETURN
      END SUBROUTINE locate_dome
!-----------------------------------------------------------------------
      SUBROUTINE set_domec
!
! ... Compute the initial conditions for a semi-spherical dome
!
      USE atmospheric_conditions, ONLY: p_atm
      USE control_flags, ONLY: job_type, lpr
      USE dimensions, ONLY: nsolid, ngas, nx
      USE domain_decomposition, ONLY: ncint, meshinds
      USE environment, ONLY: cpclock
      USE eos_gas, ONLY: ygc
      USE gas_constants, ONLY: gas_type, rgas
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_temperature, ONLY: tg, ts
      USE gas_solid_velocity, ONLY: ug, wg, vg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE grid, ONLY: flag, x, y, kv, dome_cell
      USE io_files, ONLY: testunit
      USE parallel, ONLY: mpime, root
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: ep, p
      USE array_filters, ONLY: interp
      IMPLICIT NONE

      REAL*8 :: ygcsum, ra, pi, psi, raddo, erre
      REAL*8 :: beta, p_hydro
      INTEGER :: ijk, imesh, i,j,k, is, ig, n, counter
!      
! ... Set constant for the computation of dome pressure
!
      pi = 4.D0 * ATAN(1.D0)
      psi = 2.D0 * pi
      erre = rgas / 18.D0 * 1.D3
      beta = 2.D0 * dome_gasvisc * gas_flux * erre * temperature
      beta = beta / ( p_atm(kv)**2 * permeability * psi * dome_radius )
!
      IF (lpr > 1) THEN
        IF (mpime == root) THEN
          WRITE(logunit,*) 'Dome Radial pressure profile'
          raddo = 0.D0
          DO WHILE (raddo <= dome_radius)
            WRITE(logunit,*) raddo, p_dome(raddo,p_atm(kv),beta)
            raddo = raddo + 1.D0
          END DO
        END IF
      END IF
!
! ... Set initial conditions in the dome cells
! ... (Immersed Boundaries are NOT considered within the dome!)
!
      mesh_loop: DO ijk = 1, ncint      
      IF(flag(ijk) == dome_cell) THEN
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
            ra = dcell(n)%radius
            p_hydro = dcell(n)%p_hydro
         
            ! ... Set the initial conditions, as
            ! ... specified in the input file on 
            ! ... all cells enclosing the vent
            !
            ug(ijk) = 0.D0
            IF (job_type == '3D') vg(ijk) = 0.D0 
            wg(ijk) = 0.D0
            !
            tg(ijk) = temperature
            ep(ijk) = 1.D0 - SUM(particle_fraction(1:nsolid))
            IF (idome == 1) THEN
              p(ijk)  = p_dome(ra,p_atm(kv),beta)
            ELSE IF (idome == 2) THEN
              p(ijk) = overpressure
            END IF
            !
            ! ... Add the hydrostatic pressure due to dome mass
            !
            IF (idw >= 1)  p(ijk)  = p(ijk) + p_hydro
            !
            DO ig = 1, ngas
              ygc(ijk,ig) = dome_ygc(gas_type(ig))
            END DO
            !
            DO is = 1,nsolid
              us(ijk,is)  = 0.D0
              IF (job_type == '3D') vs(ijk,is)  = 0.D0
              ws(ijk,is) = 0.D0
              !
              ts(ijk,is)  = temperature
              rlk(ijk,is) = particle_fraction(is)*rl(is)
            END DO
!          
            ! ... check gas components closure relation
            !
            ygcsum = SUM(ygc(ijk,:))
            IF ( ygcsum /= 1.D0 ) THEN
              ygc(ijk,ngas) = 1.D0 - SUM( ygc(ijk,1:ngas-1) )
            END IF
          END IF
!          
          END IF
        END DO mesh_loop
        CALL compute_dome_mass_energy

      DEALLOCATE(dcell)

      RETURN
      END SUBROUTINE set_domec
!-----------------------------------------------------------------------
      SUBROUTINE compute_dome_mass_energy
      USE atmospheric_conditions, ONLY: p_atm
      USE control_flags, ONLY: job_type, lpr
      USE dimensions, ONLY: nx, ny, nz
      USE domain_decomposition, ONLY: ncint, meshinds
      USE gas_constants, ONLY: rgas
      USE grid, ONLY: x, y, z, xb, yb, zb, r
      USE grid, ONLY: flag, dome_cell
      USE parallel, ONLY: mpime, root
      USE particles_constants, ONLY: rl
      USE pressure_epsilon, ONLY: p
      IMPLICIT NONE
      INTEGER :: ijk, i, j, k, imesh, counter, n
      REAL*8 :: cgam, gammawater, muwater
      REAL*8 :: pi, twopi
      REAL*8 :: isoe, adie, mgd, mpd, vold, dvol
      REAL*8 :: soldns, soldfr, voidfr, fact

      gammawater = 1.33D0
      muwater    = 18.D-3
      cgam = ( gammawater - 1.D0 ) / gammawater
      pi = 4.D0 * ATAN(1.D0)
      twopi = 2.D0 * pi
      fact = rgas * temperature / muwater 
      soldns = SUM( particle_fraction(1:nsolid)*rl(1:nsolid) )
      soldfr = SUM( particle_fraction(1:nsolid) )
      voidfr = 1.D0 - soldfr

      vold = 0.D0
      mgd = 0.D0
      adie = 0.D0
      isoe = 0.D0
      mesh_loop: DO ijk = 1, ncint      
        IF(flag(ijk) == dome_cell) THEN
          CALL meshinds(ijk,imesh,i,j,k)
          !
          ! ... volume of a cell
          IF (job_type == '2D') THEN
            dvol = twopi * r(i) * (xb(i)-xb(i-1))*(zb(k)-zb(k-1))
          ELSE IF (job_type == '3D') THEN
            dvol = (xb(i)-xb(i-1))*(yb(j)-yb(j-1))*(zb(k)-zb(k-1))
          END IF
          !
          ! ... specific energy 
          isoe = isoe + p(ijk) * log(p(ijk)/p_atm(k)) * dvol
          adie = adie + p(ijk) * (1.D0 - (p_atm(k)/p(ijk))**cgam) * dvol
          !
          ! ... mass per unit volume (cell)
          mgd = mgd + p(ijk) * dvol
          vold = vold + dvol
        END IF
      END DO mesh_loop

      CALL parallel_sum_real(isoe,1)
      CALL parallel_sum_real(adie,1)
      CALL parallel_sum_real(mgd,1)
      CALL parallel_sum_real(vold,1)
      !
      ! ... solid and gas mass
      mpd = soldns * vold
      mgd = mgd / fact * voidfr
      !
      ! ... 
      isoe = isoe * voidfr / (mpd + mgd)
      adie = adie / cgam * voidfr / (mpd + mgd)
!
! ... Print out the dome coordinates
!
      IF( lpr > 0 .AND. mpime == root ) THEN
        WRITE(logunit,*)
        WRITE(logunit,*) 'Mass and Energy stored in the dome'
        WRITE(logunit,100) vold 
        WRITE(logunit,200) mgd
        WRITE(logunit,201) mpd
        WRITE(logunit,300) isoe
        WRITE(logunit,400) adie
        WRITE(logunit,*)
100     FORMAT(1X,'Dome volume: ',F18.4)
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

      a1 = ACOS(r/dome_radius * COS(a))
      dh = dome_radius * SIN(a1) - r * SIN(a)

      RETURN
      END FUNCTION dh
!-----------------------------------------------------------------------
      REAL*8 FUNCTION p_dome(r2,pa,beta)
      !
      ! ... Dome pressurization model (Woods et al., 2002)
      !
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: r2, pa, beta
      REAL*8 :: r1, fact

      r1 = MAX(r2, conduit_radius)
      IF (r1 <= dome_radius) THEN
              fact = dome_radius / r1 - 1.D0
      ELSE
              fact = 0.D0
      END IF
      p_dome = pa * DSQRT(1.D0 + beta * fact)

      RETURN
      END FUNCTION p_dome
!-----------------------------------------------------------------------
      END MODULE dome_conditions
!-----------------------------------------------------------------------
