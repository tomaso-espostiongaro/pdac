!----------------------------------------------------------------------
   MODULE input_module
!----------------------------------------------------------------------

      USE dimensions, ONLY: max_nsolid, ngas, nroughx, max_size, max_nblock, max_ngas, nr, nz, nx, ny

      REAL*8 :: diameter(max_nsolid)
      REAL*8 :: density(max_nsolid)
      REAL*8 :: sphericity(max_nsolid)
      REAL*8 :: viscosity(max_nsolid)
      REAL*8 :: specific_heat(max_nsolid)
      REAL*8 :: thermal_conductivity(max_nsolid)
      REAL*8 :: dr0, dz0, dx0, dy0
      CHARACTER(LEN=80) :: run_name
      CHARACTER(LEN=80) :: restart_mode
      INTEGER :: iuni

! ... MESH
      REAL*8 :: delta_r(max_size)
      REAL*8 :: delta_x(max_size)
      REAL*8 :: delta_y(max_size)
      REAL*8 :: delta_z(max_size)

! ... FIXED_FLOWS
      INTEGER :: number_of_block
      INTEGER :: block_type(max_nblock)
      INTEGER :: block_bounds(6,max_nblock)
      REAL*8 :: fixed_vgas_r(max_nblock)
      REAL*8 :: fixed_vgas_x(max_nblock)
      REAL*8 :: fixed_vgas_y(max_nblock)
      REAL*8 :: fixed_vgas_z(max_nblock)
      REAL*8 :: fixed_pressure(max_nblock)
      REAL*8 :: fixed_gaseps(max_nblock)
      REAL*8 :: fixed_gastemp(max_nblock)
      REAL*8 :: fixed_vpart_r(max_nsolid, max_nblock)
      REAL*8 :: fixed_vpart_x(max_nsolid, max_nblock)
      REAL*8 :: fixed_vpart_y(max_nsolid, max_nblock)
      REAL*8 :: fixed_vpart_z(max_nsolid, max_nblock)
      REAL*8 :: fixed_parteps(max_nsolid, max_nblock)
      REAL*8 :: fixed_parttemp(max_nsolid, max_nblock)
      REAL*8 :: fixed_gasconc(max_ngas, max_nblock)

! ... INITIAL_CONDITIONS
      REAL*8 :: initial_vgas_r
      REAL*8 :: initial_vgas_x
      REAL*8 :: initial_vgas_y
      REAL*8 :: initial_vgas_z
      REAL*8 :: initial_pressure
      REAL*8 :: initial_void_fraction   
      REAL*8 :: max_packing          ! maximum particle packing (volumetric fraction ~0.67)
      REAL*8 :: initial_temperature
      REAL*8 :: initial_vpart_r      ! initial particle velocities (radial, x, y and z)
      REAL*8 :: initial_vpart_x
      REAL*8 :: initial_vpart_y
      REAL*8 :: initial_vpart_z
      REAL*8 :: initial_gasconc(max_ngas) ! initial gas concentration (for each specie)
!-----------------------------------------------------------------------------------------
     CONTAINS
!-----------------------------------------------------------------------------------------
!
      SUBROUTINE input( iunit )

      USE atmosphere, ONLY: w0, u0, p0, temp0, us0, ws0, ep0, epsmx0, gravx, gravz
      USE flux_limiters, ONLY: beta, muscl
      USE gas_constants, ONLY: default_gas
      USE grid, ONLY: dx, dy, dz, dr, itc, mesh_partition
      USE iterative_solver, ONLY: inmax, maxout, omega
      USE output_dump, ONLY: nfil
      USE parallel, ONLY: mpime, root
      USE particles_constants, ONLY: rl, inrl, kap, &
     &     cmus, phis, cps, dk, nsolid
      USE phases_matrix, ONLY: rlim
      USE reactions, ONLY: irex
      USE roughness_module, ONLY: zrough, allocate_roughness, roughness
      USE initial_conditions, ONLY: setup, epsob, wpob, tpob, ygc0, &
     &    ygcob, upob, wgob, ugob, pob, tgob, epob, lpr, zzero
      USE time_parameters, ONLY: time, tstop, dt, tpr, tdump, itd, & 
     &                            timestart, rungekut
      USE turbulence_model, ONLY: iturb, cmut, iss, modturbo
      USE control_flags, ONLY: job_type
      USE set_indexes, ONLY: subsc_setup
!
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: iunit

      NAMELIST / control / run_name, job_type, restart_mode, time, tstop, dt, lpr, tpr, &
        tdump, nfil, irex, iss, iturb, modturbo, cmut, rlim, gravx, gravz, &
        ngas, default_gas
!
      NAMELIST / mesh / nr, nx, ny, nz, itc, iuni, dr0, dx0, dy0, dz0, &
        mesh_partition
!
      NAMELIST / particles / nsolid, diameter, density, sphericity, &
        viscosity, specific_heat, thermal_conductivity
!
      NAMELIST / numeric / rungekut, beta, muscl, inmax, maxout, omega
!
      INTEGER :: i, j, k, n, m, ig
      CHARACTER(LEN=80) :: card
      LOGICAL :: tend
!
!    Sets default values

! ... Control

      run_name = 'run2d'
      job_type = '2D'
      restart_mode = 'from_scratch'
      time = 0.0D0
      tstop = 1.0D0
      dt = 0.01D0
      lpr = 1
      tpr = 1.0D0
      tdump = 20.0D0
      nfil = 0
      irex = 1
      iss  = 0
      iturb = 1
      modturbo = 1
      cmut = 0.1D0
      rlim = 1.0D-8
      gravx = 0.0D0
      gravz = -9.81D0
      ngas = 7
      default_gas = 6

! ... Mesh

      nr = 100
      nz = 100
      nx = 1
      ny = 1
      itc = 0
      mesh_partition = 1
      iuni = 1
      dr0  = 10.D0
      dz0  = 10.D0
      dx0  = 10.D0
      dy0  = 10.D0
      zzero = 0.0D0

! ... Particles
 
      nsolid = 2
      diameter = 100
      density = 2700
      sphericity = 1.0
      viscosity = 0.5
      specific_heat = 1.2D3
      thermal_conductivity = 2.D0

! ... Numeric

      rungekut = 1
      beta = 0.25
      muscl = 0
      inmax = 8
      maxout = 1000
      omega = 1.0

!
! reading of input file
!

      IF(mpime == root) THEN
        READ(iunit, control) 
      END IF
!
      CALL bcast_character(run_name,80,root)
      CALL bcast_character(job_type,80,root)
      CALL bcast_character(restart_mode,80,root)
      CALL bcast_real(time,1,root)
      CALL bcast_real(tstop,1,root)
      CALL bcast_real(dt,1,root)
      CALL bcast_integer(lpr,1,root)
      CALL bcast_real(tpr,1,root)
      CALL bcast_real(tdump,1,root)
      CALL bcast_integer(nfil,1,root)
      CALL bcast_integer(irex,1,root)
      CALL bcast_integer(iss,1,root)
      CALL bcast_integer(iturb,1,root)
      CALL bcast_integer(modturbo,1,root)
      CALL bcast_real(cmut,1,root)
      CALL bcast_real(rlim,1,root)
      CALL bcast_real(gravx,1,root)
      CALL bcast_real(gravz,1,root)
      CALL bcast_integer(ngas,1,root)
      CALL bcast_integer(default_gas,1,root)

      SELECT CASE ( TRIM(restart_mode) )
        CASE ('from_scratch', 'default')
          itd = 1 
        CASE ('restart')
          itd = 2
        CASE ('outp_recover')
          itd = 3
        CASE DEFAULT
          CALL error(' input ',' unknown restart_mode '//TRIM(restart_mode), 1 )
      END SELECT

      SELECT CASE ( TRIM(job_type) )
        CASE ('2D', '2d' )
          job_type = '2D'
        CASE ('3D', '3d' )
          job_type = '3D'
!          CALL error(' input ',' job_type = 3D not yet implemented ', 1 )
        CASE DEFAULT
          CALL error(' input ',' unknown job_type '//TRIM(job_type), 1 )
      END SELECT

      CALL subsc_setup( job_type )

      IF(mpime == root) THEN
        READ(iunit, mesh) 
      END IF

      CALL bcast_integer(nr,1,root)
      CALL bcast_integer(nx,1,root)
      CALL bcast_integer(ny,1,root)
      CALL bcast_integer(nz,1,root)
      CALL bcast_integer(itc,1,root)
      CALL bcast_integer(iuni,1,root)
      CALL bcast_real(dr0,1,root)
      CALL bcast_real(dx0,1,root)
      CALL bcast_real(dy0,1,root)
      CALL bcast_real(dz0,1,root)
      CALL bcast_real(zzero,1,root)

      IF(mpime == root) THEN
        READ(iunit, particles) 
      END IF

      CALL bcast_integer(nsolid,1,root)
      CALL bcast_real(diameter,max_nsolid,root)
      CALL bcast_real(density,max_nsolid,root)
      CALL bcast_real(sphericity,max_nsolid,root)
      CALL bcast_real(viscosity,max_nsolid,root)
      CALL bcast_real(specific_heat,max_nsolid,root)
      CALL bcast_real(thermal_conductivity,max_nsolid,root)

      IF( nsolid > max_nsolid .OR. nsolid < 1 ) THEN
        CALL error( ' input ', ' nsolid out of range ', nsolid )
      END IF

      IF(mpime == root) THEN
        READ(iunit, numeric) 
      END IF

      CALL bcast_integer(rungekut,1,root)
      CALL bcast_real(beta,1,root)
      CALL bcast_integer(muscl,1,root)
      CALL bcast_integer(inmax,1,root)
      CALL bcast_integer(maxout,1,root)
      CALL bcast_real(omega,1,root)

      CALL allocate_roughness( zrough, nroughx )
      tend = .FALSE.
      IF(mpime == root) THEN
        rough_search: DO 
          READ(5,*,END=100) card
          IF( TRIM(card) == 'ROUGHNESS' ) THEN
            EXIT rough_search
          END IF
        END DO rough_search
        READ(5,*) zrough%ir, (zrough%r(i),i=1,zrough%ir), zrough%roucha
        GOTO 110
 100    tend = .TRUE.
 110    continue
      END IF
!
      CALL bcast_logical(tend, 1, root)
      IF( tend ) THEN
        CALL error( ' input ', ' ROUGHNESS card not found ', 1 )
      END IF
      CALL bcast_integer(zrough%ir, 1, root)
      CALL bcast_real(zrough%r, zrough%ir, root)
      CALL bcast_real(zrough%roucha, 1, root)
!
!
      tend = .FALSE.
      IF(mpime == root) THEN
        mesh_search: DO
          READ(5,*,END=200) card
          IF( TRIM(card) == 'MESH' ) THEN
            EXIT mesh_search
          END IF
        END DO mesh_search

        IF( job_type == '2D' ) THEN
          READ(5,*) (delta_r(i),i=1,nr)
        ELSE IF( job_type == '3D' ) THEN
          READ(5,*) (delta_x(i),i=1,nx)
          READ(5,*) (delta_y(j),j=1,ny)
        END IF
        READ(5,*) (delta_z(k),k=1,nz)

        IF (iuni == 0) THEN
          delta_r = dr0
          delta_x = dx0
          delta_y = dy0
          delta_z = dz0
        ENDIF

        GOTO 210
 200    tend = .TRUE.
 210    continue
      END IF
!
      CALL bcast_logical(tend, 1, root)
      IF( tend ) THEN
        CALL error( ' input ', ' MESH card not found ', 1 )
      END IF
      CALL bcast_real(delta_r,nr,root)
      CALL bcast_real(delta_x,nx,root)
      CALL bcast_real(delta_y,ny,root)
      CALL bcast_real(delta_z,nz,root)

      tend = .FALSE.
      IF(mpime == root) THEN
        fixed_flows_search: DO
          READ(5,*,END=300) card
          IF( TRIM(card) == 'FIXED_FLOWS' ) THEN
            EXIT fixed_flows_search
          END IF
        END DO fixed_flows_search

        READ(5,*) number_of_block
        IF (job_type == '2D') THEN
          DO n = 1, number_of_block
            READ(5,*) block_type(n),block_bounds(1,n), block_bounds(2,n), block_bounds(5,n), block_bounds(6,n)
            IF( block_type(n) == 1 .OR. block_type(n) == 5) THEN
              READ(5,*) fixed_vgas_r(n), fixed_vgas_z(n), fixed_pressure(n), fixed_gaseps(n), fixed_gastemp(n)
              READ(5,*) ( fixed_vpart_r(k,n), fixed_vpart_z(k,n), fixed_parteps(k,n), fixed_parttemp(k,n), k=1, nsolid)
              READ(5,*) ( fixed_gasconc(ig,n), ig=1, ngas )
            ENDIF
          END DO
        ELSE IF (job_type == '3D') THEN
          DO n = 1, number_of_block
            READ(5,*) block_type(n), block_bounds(1,n), block_bounds(2,n), &
                                     block_bounds(3,n), block_bounds(4,n), &
                                     block_bounds(5,n), block_bounds(6,n)
            IF( block_type(n) == 1 .OR. block_type(n) == 5) THEN
              READ(5,*) fixed_vgas_x(n), fixed_vgas_y(n), fixed_vgas_z(n), &
                        fixed_pressure(n), fixed_gaseps(n), fixed_gastemp(n)
              READ(5,*) ( fixed_vpart_x(k,n), fixed_vpart_y(k,n), fixed_vpart_z(k,n), &
                          fixed_parteps(k,n), fixed_parttemp(k,n), k=1, nsolid)
              READ(5,*) ( fixed_gasconc(ig,n), ig=1, ngas )
            ENDIF
          END DO
        ELSE
          CALL error( ' input # FIXED FLOW ', ' unknown job_type', 1 )
        END IF

        GOTO 310
 300    tend = .TRUE.
 310    continue
      END IF
!
      CALL bcast_logical(tend, 1, root)
      IF( tend ) THEN
        CALL error( ' input ', ' FIXED FLOWS card not found ', 1 )
      END IF
      CALL bcast_integer(number_of_block, 1, root)
      CALL bcast_integer(block_type, SIZE(block_type), root)
      CALL bcast_integer(block_bounds, SIZE(block_bounds), root)
      CALL bcast_real(fixed_vgas_r, SIZE(fixed_vgas_r), root)
      CALL bcast_real(fixed_vgas_x, SIZE(fixed_vgas_x), root)
      CALL bcast_real(fixed_vgas_y, SIZE(fixed_vgas_y), root)
      CALL bcast_real(fixed_vgas_z, SIZE(fixed_vgas_z), root)
      CALL bcast_real(fixed_pressure, SIZE(fixed_pressure), root)
      CALL bcast_real(fixed_gaseps, SIZE(fixed_gaseps), root)
      CALL bcast_real(fixed_gastemp, SIZE(fixed_gastemp), root)
      CALL bcast_real(fixed_vpart_r, SIZE(fixed_vpart_r), root)
      CALL bcast_real(fixed_vpart_x, SIZE(fixed_vpart_x), root)
      CALL bcast_real(fixed_vpart_y, SIZE(fixed_vpart_y), root)
      CALL bcast_real(fixed_vpart_z, SIZE(fixed_vpart_z), root)
      CALL bcast_real(fixed_parteps, SIZE(fixed_parteps), root)
      CALL bcast_real(fixed_parttemp, SIZE(fixed_parttemp), root)
      CALL bcast_real(fixed_gasconc, SIZE(fixed_gasconc), root)


      tend = .FALSE.
      IF(mpime == root) THEN
        initial_conditions_search: DO
          READ(5,*,END=400) card
          IF( TRIM(card) == 'INITIAL_CONDITIONS' ) THEN
            EXIT initial_conditions_search
          END IF
        END DO initial_conditions_search

        IF( job_type == '2D' ) THEN
          READ(5,*) initial_vgas_r, initial_vgas_z, initial_pressure, initial_void_fraction, max_packing, &
            initial_temperature
          READ(5,*) initial_vpart_r, initial_vpart_z
        ELSE IF( job_type == '3D' ) THEN
          READ(5,*) initial_vgas_x, initial_vgas_y, initial_vgas_z, initial_pressure, initial_void_fraction, &
            max_packing, initial_temperature
          READ(5,*) initial_vpart_x, initial_vpart_y, initial_vpart_z
        ELSE 
          CALL error('input # INITIAL_CONDITIONS', 'unknown job_type',1) 
        ENDIF
        READ(5,*) (initial_gasconc(ig), ig=1, ngas)

        GOTO 410
 400    tend = .TRUE.
 410    continue
      END IF
!
      CALL bcast_logical(tend, 1, root)
      IF( tend ) THEN
        CALL error( ' input ', ' INITIAL_CONDITIONS card not found ', 1 )
      END IF
      CALL bcast_real(initial_vgas_r,1,root)
      CALL bcast_real(initial_vgas_x,1,root)
      CALL bcast_real(initial_vgas_y,1,root)
      CALL bcast_real(initial_vgas_z,1,root)
      CALL bcast_real(initial_pressure,1,root)
      CALL bcast_real(initial_void_fraction,1,root)
      CALL bcast_real(max_packing,1,root)
      CALL bcast_real(initial_temperature,1,root)
      CALL bcast_real(initial_vpart_r,1,root)
      CALL bcast_real(initial_vpart_x,1,root)
      CALL bcast_real(initial_vpart_y,1,root)
      CALL bcast_real(initial_vpart_z,1,root)
      CALL bcast_real(initial_gasconc, SIZE(initial_gasconc),root)

      END SUBROUTINE

!----------------------------------------------------------------------
      END MODULE input_module
!----------------------------------------------------------------------
