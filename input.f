!----------------------------------------------------------------------
   MODULE input_module
!----------------------------------------------------------------------

      USE dimensions, ONLY: max_nsolid, ngas, nroughx, max_size, max_nblock, max_ngas, nz, nx, ny

      REAL*8 :: diameter(max_nsolid)
      REAL*8 :: density(max_nsolid)
      REAL*8 :: sphericity(max_nsolid)
      REAL*8 :: viscosity(max_nsolid)
      REAL*8 :: specific_heat(max_nsolid)
      REAL*8 :: thermal_conductivity(max_nsolid)
      REAL*8 :: dz0, dx0, dy0
      CHARACTER(LEN=80) :: run_name
      CHARACTER(LEN=80) :: restart_mode
      INTEGER :: iuni

! ... MESH
      REAL*8 :: delta_x(max_size)
      REAL*8 :: delta_y(max_size)
      REAL*8 :: delta_z(max_size)

      REAL*8 :: origin_x
      REAL*8 :: origin_y
      REAL*8 :: origin_z

! ... FIXED_FLOWS
      INTEGER :: number_of_block
      INTEGER :: block_type(max_nblock)
      INTEGER :: block_bounds(6,max_nblock)
      REAL*8 :: fixed_vgas_x(max_nblock)
      REAL*8 :: fixed_vgas_y(max_nblock)
      REAL*8 :: fixed_vgas_z(max_nblock)
      REAL*8 :: fixed_pressure(max_nblock)
      REAL*8 :: fixed_gaseps(max_nblock)
      REAL*8 :: fixed_gastemp(max_nblock)
      REAL*8 :: fixed_vpart_x(max_nsolid, max_nblock)
      REAL*8 :: fixed_vpart_y(max_nsolid, max_nblock)
      REAL*8 :: fixed_vpart_z(max_nsolid, max_nblock)
      REAL*8 :: fixed_parteps(max_nsolid, max_nblock)
      REAL*8 :: fixed_parttemp(max_nsolid, max_nblock)
      REAL*8 :: fixed_gasconc(max_ngas, max_nblock)

! ... INITIAL_CONDITIONS
      REAL*8 :: initial_vgas_x
      REAL*8 :: initial_vgas_y
      REAL*8 :: initial_vgas_z
      REAL*8 :: initial_pressure
      REAL*8 :: initial_void_fraction   
      REAL*8 :: max_packing          ! maximum particle packing (volumetric fraction ~0.67)
      REAL*8 :: initial_temperature
      REAL*8 :: initial_vpart_x      ! initial particle velocities (x, y and z)
      REAL*8 :: initial_vpart_y
      REAL*8 :: initial_vpart_z
      REAL*8 :: initial_gasconc(max_ngas) ! initial gas concentration (for each specie)

      INTEGER :: iuni_nml = 20
      INTEGER :: iuni_fld = 21
      INTEGER :: iuni_grx = 22
      INTEGER :: iuni_gry = 23
      INTEGER :: iuni_grz = 24

! ... POST PROCESSING
      INTEGER :: first_out = 1
      INTEGER :: last_out = 1
      INTEGER :: incr_out = 1

!-----------------------------------------------------------------------
     CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE input( iunit, which )

      USE atmosphere, ONLY: w0, u0, p0, temp0, us0, ws0, ep0, epsmx0, gravx, gravz
      USE control_flags, ONLY: nfil, job_type
      USE domain_decomposition, ONLY: mesh_partition
      USE flux_limiters, ONLY: beta, muscl, lim_type
      USE gas_constants, ONLY: default_gas
      USE gas_solid_viscosity, ONLY: gas_viscosity, part_viscosity
      USE grid, ONLY: dx, dy, dz, itc
      USE initial_conditions, ONLY: setup, epsob, wpob, tpob, ygc0, &
     &    ygcob, upob, wgob, ugob, pob, tgob, epob, lpr, zzero
      USE iterative_solver, ONLY: inmax, maxout, omega
      USE io_restart, ONLY: old_restart, max_seconds
      USE output_dump, ONLY: formatted_output

      USE parallel, ONLY: mpime, root
      USE particles_constants, ONLY: rl, inrl, kap, &
     &     cmus, phis, cps, dk, nsolid
      USE phases_matrix, ONLY: rlim
      USE reactions, ONLY: irex
      USE roughness_module, ONLY: zrough, allocate_roughness, roughness
      USE set_indexes, ONLY: subsc_setup
      USE time_parameters, ONLY: time, tstop, dt, tpr, tdump, itd, & 
     &                            timestart, rungekut
      USE turbulence_model, ONLY: iturb, cmut, iss, modturbo
!
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: iunit
      CHARACTER(LEN=2) :: which

      NAMELIST / control / run_name, job_type, restart_mode,       &
        time, tstop, dt, lpr, tpr, tdump, nfil, irex,              &
        gas_viscosity, part_viscosity, iss, iturb, modturbo, cmut, &
        rlim, gravx, gravz, ngas, default_gas,                     &
        formatted_output, old_restart, max_seconds

      NAMELIST / pp / first_out, last_out, incr_out
!
      NAMELIST / mesh / nx, ny, nz, itc, iuni, dx0, dy0, dz0, &
        origin_x, origin_y, origin_z, mesh_partition
!
      NAMELIST / particles / nsolid, diameter, density, sphericity, &
        viscosity, specific_heat, thermal_conductivity
!
      NAMELIST / numeric / rungekut, beta, muscl, lim_type, inmax, maxout, omega
!
      INTEGER :: i, j, k, n, m, ig
      REAL*8 :: grx, gry, grz
      CHARACTER(LEN=80) :: card
      LOGICAL :: tend
!
!    Sets default values

! ... Control

      run_name = 'run2d'
      job_type = '2D'
      restart_mode = 'from_scratch'
      time = 0.0D0      ! start time in seconds
      tstop = 1.0D0     ! stop time seconds
      dt = 0.01D0       ! time increment seconds
      lpr = 1           ! verbosity (not yet implemented)
      tpr = 1.0D0       ! write to output file every tpr seconds of simulated time
      tdump = 20.0D0    ! write restart every tdump seconds of simulated time
      nfil = 0          ! output file index
      irex = 1          ! ( 1 no reaction, 2 use hrex. Not used )
      gas_viscosity  = .TRUE. ! include molecular gas viscosity
      part_viscosity = .TRUE. ! include collisional particle viscosity
      iss  = 0          ! ( 0 no solid turbulent stress, 1 compute stress )
      iturb = 1         ! turbulence  ( 0 no turbo, 1 turbo, 2 turbo + rough )
      modturbo = 1      ! turbulence  ( 1 smag, 2 dynamic )
      cmut = 0.1D0      ! Smagorinsky constant
      rlim = 1.0D-8     ! limit for off-diagonal contribution in matrix inversion
      gravx = 0.0D0     ! gravity along x
      gravz = -9.81D0   ! gravity along z
      ngas = 7          ! number of gas phase
      default_gas = 6   ! atmosphere 
      formatted_output = .TRUE.
      old_restart = .FALSE.
      max_seconds = 20000.0

! ... Mesh

      nx = 100                !  number of cell in the X directions
      nz = 100                !  number of cell in the Z directions
      ny = 1                  !  number of cell in the Y directions
      itc = 0                 !  itc = 1 cylindrical coordinates are used
      mesh_partition = 1      !  type of partition ( 1 = layers, 2 = columns, 3 = blocks )
      iuni = 0                !  1 = uniform grid, 0 = non uniform grid
      dz0  = 10.D0            !  default cell z size in meters
      dx0  = 10.D0            !  default cell x size in meters
      dy0  = 10.D0            !  default cell y size in meters
      origin_x  = 0.D0        !  default x coo. of the system origin
      origin_y  = 0.D0        !  default y coo. of the system origin
      origin_z  = 0.D0        !  default z coo. of the system origin
      zzero = 0.0D0           !  grid bottom level

! ... Particles
 
      nsolid = 2              !  number of solid components
      diameter = 100          !  particle diameters in micron
      density = 2700          !  particle density in kg/m^3
      sphericity = 1.0        !  sphericity coefficients ( 1.0 sphere )
      viscosity = 0.5         !  viscosity coefficient ( Pa * sec. )
      specific_heat = 1.2D3   !  specific heat ( Joule / ( Kelvin * Kg ) )
      thermal_conductivity = 2.D0 ! thermal_conductivity ( ... ) 

! ... Numeric

      rungekut = 1      !  number of runge-kutta cycles
      beta = 0.25       !  upwinding degree
      lim_type = 0      !  limiter type 
      muscl = 0         !  0 first order, 1 muscl ( high order )
      inmax = 8         !  maximum number of pressure correction steps
      maxout = 1000     !  maximum number of solver iteration
      omega = 1.0       !  relaxation parameter  ( 0.5 under - 2.0 over)

!
! reading of input file
!

      IF( which == 'PP' ) THEN 
        IF( mpime == root ) THEN
          OPEN( UNIT=iuni_nml, FILE='pdac.nml', STATUS='UNKNOWN')
          OPEN( UNIT=iuni_grx, FILE='pdac.grx', STATUS='UNKNOWN')
          OPEN( UNIT=iuni_gry, FILE='pdac.gry', STATUS='UNKNOWN')
          OPEN( UNIT=iuni_grz, FILE='pdac.grz', STATUS='UNKNOWN')
          OPEN( UNIT=iuni_fld, FILE='pdac.fld', STATUS='UNKNOWN')
        END IF
      END IF

      IF(mpime == root) THEN
        READ(iunit, control) 
        IF( which == 'PP' ) THEN
          WRITE(iuni_nml, control) 
        END IF
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
      CALL bcast_logical(formatted_output,1,root)
      CALL bcast_real(max_seconds,1,root)
      CALL bcast_logical(gas_viscosity,1,root)
      CALL bcast_logical(part_viscosity,1,root)

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
        IF( which == 'PP' ) THEN
          WRITE(iuni_nml, mesh) 
        END IF
      END IF

      CALL bcast_integer(nx,1,root)
      CALL bcast_integer(ny,1,root)
      CALL bcast_integer(nz,1,root)
      CALL bcast_integer(itc,1,root)
      CALL bcast_integer(iuni,1,root)
      CALL bcast_real(dx0,1,root)
      CALL bcast_real(dy0,1,root)
      CALL bcast_real(dz0,1,root)
      CALL bcast_real(origin_x,1,root)
      CALL bcast_real(origin_y,1,root)
      CALL bcast_real(origin_z,1,root)
      CALL bcast_real(zzero,1,root)

      IF(mpime == root) THEN
        READ(iunit, particles) 
        IF( which == 'PP' ) THEN
          WRITE(iuni_nml, particles) 
        END IF
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
        IF( which == 'PP' ) THEN
          WRITE(iuni_nml, numeric) 
        END IF
      END IF

      CALL bcast_integer(rungekut,1,root)
      CALL bcast_real(beta,1,root)
      CALL bcast_integer(lim_type,1,root)
      CALL bcast_integer(muscl,1,root)
      CALL bcast_integer(inmax,1,root)
      CALL bcast_integer(maxout,1,root)
      CALL bcast_real(omega,1,root)

      IF( which == 'PP' ) THEN
        IF(mpime == root) THEN
          READ(iunit, pp) 
          WRITE(iuni_nml, pp) 
        END IF
        CALL bcast_integer(first_out,1,root)
        CALL bcast_integer(last_out,1,root)
        CALL bcast_integer(incr_out,1,root)
      END IF

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
!     Read roughness parameters
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
!     Read cell sizes
!
      tend = .FALSE.
      IF(mpime == root) THEN
        mesh_search: DO
          READ(5,*,END=200) card
          IF( TRIM(card) == 'MESH' ) THEN
            EXIT mesh_search
          END IF
        END DO mesh_search

        IF (iuni == 0) THEN
          READ(5,*) (delta_x(i),i=1,nx)
          IF( job_type == '3D' ) THEN
            READ(5,*) (delta_y(j),j=1,ny)
          END IF
          READ(5,*) (delta_z(k),k=1,nz)
        ELSE IF (iuni == 1) THEN
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
      CALL bcast_real(delta_x,nx,root)
      CALL bcast_real(delta_y,ny,root)
      CALL bcast_real(delta_z,nz,root)

      IF(mpime == root) THEN

        IF( which == 'PP' ) THEN

          grx = origin_x - delta_x(1)
          grx = grx + delta_x(1)
          DO i = 2, nx
            WRITE( iuni_grx, * ) grx
            grx = grx + delta_x(i)
          END DO
          ! WRITE( iuni_grx, * ) grx
  
          gry = origin_y - delta_y(1)
          gry = gry + delta_y(1)
          DO i = 2, ny
            WRITE( iuni_gry, * ) gry
            gry = gry + delta_y(i)
          END DO
          ! WRITE( iuni_gry, * ) gry
  
          grz = origin_z - delta_z(1)
          grz = grz + delta_z(1)
          DO i = 2, nz
            WRITE( iuni_grz, * ) grz
            grz = grz + delta_z(i)
          END DO
          ! WRITE( iuni_grz, * ) grz

          !WRITE( iuni_fld, fmt = "('# AVS field file')" )
          !WRITE( iuni_fld, fmt = "('ndim=',I3)" ) 3
          !WRITE( iuni_fld, fmt = "('dim1=',I3)" ) nx+1
          !WRITE( iuni_fld, fmt = "('dim2=',I3)" ) ny+1
          !WRITE( iuni_fld, fmt = "('dim3=',I3)" ) nz+1
          !WRITE( iuni_fld, fmt = "('nspace=',I3)" ) 3
          !WRITE( iuni_fld, fmt = "('veclen=',I3)" ) 3
          !WRITE( iuni_fld, fmt = "('data=double')" )
          !WRITE( iuni_fld, fmt = "('field=irregular')" )
          !WRITE( iuni_fld, fmt = "('coord 1 file=pdac.grx, filetype=ascii')" )
          !WRITE( iuni_fld, fmt = "('coord 2 file=pdac.gry, filetype=ascii')" )
          !WRITE( iuni_fld, fmt = "('coord 3 file=pdac.grz, filetype=ascii')" )

        END IF

      END IF

      

!
!     Read boundary conditions
!

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
              READ(5,*) fixed_vgas_x(n), fixed_vgas_z(n), fixed_pressure(n), fixed_gaseps(n), fixed_gastemp(n)
              READ(5,*) ( fixed_vpart_x(k,n), fixed_vpart_z(k,n), fixed_parteps(k,n), fixed_parttemp(k,n), k=1, nsolid)
              READ(5,*) ( fixed_gasconc(ig,n), ig=1, ngas )
            ELSE IF( block_type(n) == 7) THEN
              OPEN(UNIT=17, FILE='input.prof', STATUS='UNKNOWN')
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
      CALL bcast_real(fixed_vgas_x, SIZE(fixed_vgas_x), root)
      CALL bcast_real(fixed_vgas_y, SIZE(fixed_vgas_y), root)
      CALL bcast_real(fixed_vgas_z, SIZE(fixed_vgas_z), root)
      CALL bcast_real(fixed_pressure, SIZE(fixed_pressure), root)
      CALL bcast_real(fixed_gaseps, SIZE(fixed_gaseps), root)
      CALL bcast_real(fixed_gastemp, SIZE(fixed_gastemp), root)
      CALL bcast_real(fixed_vpart_x, SIZE(fixed_vpart_x), root)
      CALL bcast_real(fixed_vpart_y, SIZE(fixed_vpart_y), root)
      CALL bcast_real(fixed_vpart_z, SIZE(fixed_vpart_z), root)
      CALL bcast_real(fixed_parteps, SIZE(fixed_parteps), root)
      CALL bcast_real(fixed_parttemp, SIZE(fixed_parttemp), root)
      CALL bcast_real(fixed_gasconc, SIZE(fixed_gasconc), root)

!
!     Read initial conditions ( atmosphere conditions )
!

      tend = .FALSE.
      IF(mpime == root) THEN
        initial_conditions_search: DO
          READ(5,*,END=400) card
          IF( TRIM(card) == 'INITIAL_CONDITIONS' ) THEN
            EXIT initial_conditions_search
          END IF
        END DO initial_conditions_search

        IF( job_type == '2D' ) THEN
          READ(5,*) initial_vgas_x, initial_vgas_z, initial_pressure, initial_void_fraction, max_packing, &
            initial_temperature
          READ(5,*) initial_vpart_x, initial_vpart_z
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
      CALL bcast_real(initial_vgas_x,1,root)
      CALL bcast_real(initial_vgas_y,1,root)
      CALL bcast_real(initial_vgas_z,1,root)
      CALL bcast_real(initial_pressure,1,root)
      CALL bcast_real(initial_void_fraction,1,root)
      CALL bcast_real(max_packing,1,root)
      CALL bcast_real(initial_temperature,1,root)
      CALL bcast_real(initial_vpart_x,1,root)
      CALL bcast_real(initial_vpart_y,1,root)
      CALL bcast_real(initial_vpart_z,1,root)
      CALL bcast_real(initial_gasconc, SIZE(initial_gasconc),root)

      IF( mpime == root ) THEN
        IF( which == 'PP' ) THEN
          CLOSE( UNIT=iuni_nml )
          CLOSE( UNIT=iuni_fld )
          CLOSE( UNIT=iuni_grx )
          CLOSE( UNIT=iuni_gry )
          CLOSE( UNIT=iuni_grz )
        END IF
      END IF

      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE initc

      USE atmosphere, ONLY: v0, u0, w0, p0, temp0, us0, vs0, ws0, &
     &                      ep0, epsmx0, gravx, gravy, gravz
      USE dimensions
      USE gas_constants, ONLY: present_gas
      USE grid, ONLY: dx, dy, dz, itc
      USE grid, ONLY: iob
      USE initial_conditions, ONLY: epsob, tpob, ygc0, ygcob,   &
     &     ugob, vgob, wgob, upob, vpob, wpob, pob, tgob, epob, &
     &     lpr, zzero
      USE particles_constants, ONLY: rl, inrl, kap, cmus, phis, cps, dk
      USE control_flags, ONLY: job_type

      dx(1:nx) = delta_x(1:nx)
      IF( job_type == '3D' ) THEN
        dy(1:ny) = delta_y(1:ny)
      ELSE IF( job_type == '2D' ) THEN
        dy(1:ny) = 0.D0
      END IF
      dz(1:nz) = delta_z(1:nz)
!
      iob(1:no)%typ = block_type(1:no)

      iob(1:no)%xlo = block_bounds(1,1:no)
      iob(1:no)%xhi = block_bounds(2,1:no)
      IF ( job_type == '3D' ) THEN
        iob(1:no)%ylo = block_bounds(3,1:no)
        iob(1:no)%yhi = block_bounds(4,1:no)
      END IF
      iob(1:no)%zlo = block_bounds(5,1:no)
      iob(1:no)%zhi = block_bounds(6,1:no)

      ugob(1:no)  = fixed_vgas_x(1:no)
      IF( job_type == '3D' ) THEN
        vgob(1:no)  = fixed_vgas_y(1:no)
      END IF
      wgob(1:no)  = fixed_vgas_z(1:no)
      pob(1:no)  = fixed_pressure(1:no)
      epob(1:no)  = fixed_gaseps(1:no)
      tgob(1:no)  = fixed_gastemp(1:no)
      upob(1:nsolid,1:no) = fixed_vpart_x(1:nsolid,1:no)
      IF( job_type == '3D' ) THEN
        vpob(1:nsolid,1:no) = fixed_vpart_y(1:nsolid,1:no)
      END IF
      wpob(1:nsolid,1:no) = fixed_vpart_z(1:nsolid,1:no)
      epsob(1:nsolid,1:no) = fixed_parteps(1:nsolid,1:no)
      tpob(1:nsolid,1:no) = fixed_parttemp(1:nsolid,1:no)
      ygcob(1:ngas,1:no) = fixed_gasconc(1:ngas,1:no)

      u0 = initial_vgas_x
      IF( job_type == '3D' ) THEN
        v0 = initial_vgas_y
      END IF
      w0 = initial_vgas_z

      p0 = initial_pressure
      ep0 = initial_void_fraction
      epsmx0 = max_packing
      temp0 = initial_temperature

      us0 = initial_vpart_x
      IF( job_type == '3D' ) THEN
        vs0 = initial_vpart_y
      END IF
      ws0 = initial_vpart_z

      ygc0(1:ngas) = initial_gasconc(1:ngas)

      DO ig = 1, ngas
        present_gas(ig) = (ygc0(ig) /= 0.D0 .OR. ANY(ygcob(ig,:) /= 0.D0 ) )
      END DO
!
      dk(1:nsolid) = diameter(1:nsolid) * 1.D-6 ! input diameter in microns !
      rl(1:nsolid) = density(1:nsolid)
      phis(1:nsolid) = sphericity(1:nsolid)
      cmus(1:nsolid) = viscosity(1:nsolid)
      cps(1:nsolid) = specific_heat(1:nsolid)
      kap(1:nsolid) = thermal_conductivity(1:nsolid)
!
      inrl(:)=1.D0/rl(:)

      END SUBROUTINE initc
!----------------------------------------------------------------------
      END MODULE input_module
!----------------------------------------------------------------------
