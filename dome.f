!-----------------------------------------------------------------------
      MODULE dome_conditions
!-----------------------------------------------------------------------
      USE dimensions, ONLY: max_nsolid, max_ngas, nsolid, ngas
      USE io_files, ONLY: errorunit, logunit
      IMPLICIT NONE

      ! ... flags
      !
      INTEGER :: idome

      REAL*8 :: xdome, ydome, dome_volume, dome_radius
      REAL*8 :: temperature
      REAL*8 :: particle_fraction(max_nsolid)
      !
      REAL*8 :: dome_ygc(max_ngas)

      TYPE dome_cell
        INTEGER :: imesh
        REAL*8  :: radius
      END TYPE dome_cell

      TYPE(dome_cell), ALLOCATABLE :: dcell(:)

      INTEGER :: ndm
      REAL*8 :: gas_flux
      REAL*8 :: permeability
      REAL*8 :: dome_gasvisc
      PUBLIC
      
      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
! ... This routine locate the vent on a georeferenced mesh and set the
! ... cell flags onto the neighbours of the inlet
!
      SUBROUTINE locate_dome

      USE control_flags, ONLY: job_type, lpr
      USE dimensions, ONLY: nx, ny, nz, ntot
      USE grid, ONLY: x, y, z, fl, xb, yb, zb, dz
      USE grid, ONLY: bottom, iv, jv, kv, grigen, flag
      USE volcano_topography, ONLY: itp, ord2d
      USE volcano_topography, ONLY: flatten_dem
      USE parallel, ONLY: mpime, root

      IMPLICIT NONE
      
      INTEGER :: i, j, k, ijk, n
      REAL*8 :: distance2, dvol, pi
      
      IF( job_type == '2D') RETURN
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
! ... Compute the dome radius, given the volume.
! ... The dome shape is half-a-sphere. 
!
      pi = 4.D0 * ATAN(1.D0)
      dome_radius = ( 1.5D0 * dome_volume / pi )**(1.0/3.0)
!
! ... Define the quota of the volcanic dome
! ... (considering the topography). If the topography
! ... has been red from a file, 'kv' has already been set
! ... in the 'volcano_topography' module.
!
      IF( itp < 1 ) THEN
        DO k= 1, nz
          ijk = iv + (jv-1) * nx + (k-1) * nx * ny
          IF (fl(ijk) == 3) kv = k
        END DO
      END IF
!
! ... Count the cells of the dome
!
      ndm = 0
      DO k = 1, nz
        DO j = 1, ny
          DO i = 1, nx
            ijk = i + (j-1) * nx + (k-1) * nx * ny
            distance2 = (x(i)-x(iv))**2 + (y(j)-y(jv))**2 + (z(k)-z(kv))**2
            IF (distance2 <= dome_radius**2 .AND. fl(ijk) == 1 ) THEN
                    ndm = ndm + 1
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
            distance2 = (x(i)-x(iv))**2 + (y(j)-y(jv))**2 + (z(k)-z(kv))**2
            IF (distance2 <= dome_radius**2 .AND. fl(ijk) == 1 ) THEN
                    n = n + 1
                    dcell(n)%imesh = ijk
                    dcell(n)%radius = DSQRT(distance2)
            END IF
          END DO
        END DO
      END DO
!      
! ... Compute the volume of the dome
!
      dvol = 0.D0
      DO n = 1, ndm
        ijk = dcell(n)%imesh
        k = ( ijk - 1 ) / ( nx*ny ) + 1
        j = MOD( ijk - 1, nx*ny) / nx + 1
        i = MOD( MOD( ijk - 1, nx*ny ), nx ) + 1
        dvol = dvol + (xb(i)-xb(i-1))*(yb(j)-yb(j-1))*(zb(k)-zb(k-1))
      END DO
!
! ... Print out the dome coordinates and volume
!
      IF( lpr > 0 .AND. mpime == root ) THEN
        WRITE(logunit,100) iv, jv, kv
        WRITE(logunit,200) x(iv), y(jv), z(kv)
        WRITE(logunit,300) dvol
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
      USE dimensions, ONLY: nsolid, ngas
      USE domain_decomposition, ONLY: ncint, meshinds
      USE environment, ONLY: cpclock
      USE eos_gas, ONLY: ygc
      USE gas_constants, ONLY: gas_type, rgas
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_temperature, ONLY: tg, ts
      USE gas_solid_velocity, ONLY: ug, wg, vg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE grid, ONLY: flag, x, y, kv
      USE parallel, ONLY: mpime, root
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: ep, p
      USE array_filters, ONLY: interp
      IMPLICIT NONE

      REAL*8 :: ygcsum, ra, pi, psi, raddo, erre
      REAL*8 :: beta
      INTEGER :: ijk, imesh, i,j,k, is, ig, n, counter

      IF (job_type == '2D') RETURN
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

      mesh_loop: DO ijk = 1, ncint      
        IF(flag(ijk) == 1) THEN
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
            p(ijk)  = p_dome(ra,p_atm(kv),beta)
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

      DEALLOCATE(dcell)

      RETURN
      END SUBROUTINE set_domec
!-----------------------------------------------------------------------
      REAL*8 FUNCTION p_dome(r2,pa,beta)
      !
      ! ... Dome pressurization model (Woods et al., 2002)
      !
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: r2, pa, beta
      REAL*8 :: r1, fact

      r1 = MAX(r2, 10.D0)
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
