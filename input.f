!----------------------------------------------------------------------
   MODULE input_module
!----------------------------------------------------------------------

      USE dimensions, ONLY: max_nsolid, ngas, nroughx, max_size, &
          max_nblock, max_ngas, nz, nx, ny
      USE iotk_module

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

! ... SPECIFIED_PROFILE
      REAL*8 :: prof_vgas_x(max_size)
      REAL*8 :: prof_vgas_y(max_size)
      REAL*8 :: prof_vgas_z(max_size)
      REAL*8 :: prof_pressure(max_size)
      REAL*8 :: prof_gaseps(max_size)
      REAL*8 :: prof_gastemp(max_size)
      REAL*8 :: prof_vpart_x(max_nsolid, max_size)
      REAL*8 :: prof_vpart_y(max_nsolid, max_size)
      REAL*8 :: prof_vpart_z(max_nsolid, max_size)
      REAL*8 :: prof_parteps(max_nsolid, max_size)
      REAL*8 :: prof_parttemp(max_nsolid, max_size)
      REAL*8 :: prof_gasconc(max_ngas, max_size)

! ... INITIAL_CONDITIONS
      REAL*8 :: initial_vgas_x
      REAL*8 :: initial_vgas_y
      REAL*8 :: initial_vgas_z
      REAL*8 :: initial_pressure
      REAL*8 :: initial_void_fraction   
      REAL*8 :: max_packing
      REAL*8 :: initial_temperature
      REAL*8 :: initial_vpart_x
      REAL*8 :: initial_vpart_y
      REAL*8 :: initial_vpart_z
      REAL*8 :: initial_gasconc(max_ngas)

      INTEGER :: iuni_nml = 20
      INTEGER :: iuni_fld = 21
      INTEGER :: iuni_grx = 22
      INTEGER :: iuni_gry = 23
      INTEGER :: iuni_grz = 24

! ... POST PROCESSING
      INTEGER :: first_out = 1
      INTEGER :: last_out = 1
      INTEGER :: incr_out = 1

      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE input( iunit, which )

      USE atmosphere, ONLY: w0, u0, p0, temp0, us0, ws0, ep0, epsmx0
      USE atmosphere, ONLY: gravx, gravy, gravz
      USE control_flags, ONLY: nfil, job_type, lpr, immb, itp
      USE control_flags, ONLY: implicit_fluxes, implicit_enthalpy
      USE domain_decomposition, ONLY: mesh_partition
      USE flux_limiters, ONLY: beta, muscl, lim_type
      USE gas_solid_viscosity, ONLY: gas_viscosity, part_viscosity
      USE gas_solid_viscosity, ONLY: repulsive_model
      USE grid, ONLY: dx, dy, dz, itc, zzero
      USE grid, ONLY: west, east, south, north, bottom, top, topography
      USE initial_conditions, ONLY: setup, epsob, wpob, tpob, ygc0, &
     &    ygcob, upob, wgob, ugob, pob, tgob, epob, density_specified
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
      USE initial_conditions, ONLY: ugpr, wgpr, ppr, eppr, tgpr, &
                                   uppr, wppr, epspr, tppr, ygcpr, npr
!
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: iunit
      CHARACTER(LEN=2) :: which

      NAMELIST / control / run_name, job_type, restart_mode,       &
        time, tstop, dt, lpr, tpr, tdump, nfil,                    &
        formatted_output, old_restart, max_seconds

      NAMELIST / model / irex, gas_viscosity, part_viscosity,      &
        iss, repulsive_model, iturb, modturbo, cmut, rlim,         &
        gravx, gravy, gravz, ngas, density_specified

      NAMELIST / pp / first_out, last_out, incr_out

      NAMELIST / mesh / nx, ny, nz, itc, iuni, dx0, dy0, dz0, &
        origin_x, origin_y, origin_z, mesh_partition

      NAMELIST / boundaries / west, east, south, north, bottom, top, &
        itp, topography, immb

      NAMELIST / particles / nsolid, diameter, density, sphericity, &
        viscosity, specific_heat, thermal_conductivity

      NAMELIST / numeric / rungekut, beta, muscl, lim_type, &
        inmax, maxout, omega, implicit_fluxes, implicit_enthalpy

      INTEGER :: i, j, k, n, m, ig, ierr
      REAL*8, ALLOCATABLE :: grx(:), gry(:), grz(:)
      CHARACTER(LEN=80) :: card
      CHARACTER(LEN=256) :: attr
      LOGICAL :: tend
!
!:::::::::::::::::::::::::::::  Sets default values  ::::::::::::::::::::::::::
!

! ... Control

      run_name = 'pdac_run_2d'
      job_type = '2D'
      restart_mode = 'from_scratch'
      time = 0.0D0      ! start time in seconds
      tstop = 100.0D0     ! stop time seconds
      dt = 0.01D0       ! time increment seconds
      lpr = 2           ! verbosity (not yet implemented)
      tpr = 1.0D0       ! write to output file every tpr seconds of simulated time
      tdump = 20.0D0    ! write restart every tdump seconds of simulated time
      nfil = 0          ! output file index
      formatted_output = .TRUE.
      old_restart = .FALSE.
      max_seconds = 20000.0

! ... Model

      irex = 1          ! ( 1 no reaction, 2 use hrex. Not used )
      gas_viscosity  = .TRUE. ! include molecular gas viscosity
      part_viscosity = .TRUE. ! include collisional particle viscosity
      density_specified = .FALSE. ! density specified instead of temperature
      iss  = 0          ! ( 0 no solid turbulent stress, 1 sub-grid stress) 
      repulsive_model = 1 ! ( 0 no Coulombic repulsive model )
      iturb = 1         ! turbulence  ( 0 no turbo, 1 turbo, 2 turbo + rough )
      modturbo = 1      ! turbulence  ( 1 smag, 2 dynamic )
      cmut = 0.1D0      ! Smagorinsky constant
      rlim = 1.0D-10     ! limit for off-diagonal contribution in matrix 
                        ! inversion
      gravx = 0.0D0     ! gravity along x
      gravy = 0.0D0     ! gravity along y
      gravz = -9.81D0   ! gravity along z
      ngas = 2          ! max number of gas components

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

! ... Boundaries

      west = 2                !
      east = 4                !
      bottom = 3              ! boundary types
      top = 4                 !
      south = 2               !
      north = 2               !
      itp   = 0               ! itp = 1 => read topography from file
      topography = 'topo.dat' ! file containing the topographic profile
      immb  = 0               ! 1: use immersed boundaries

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
      lim_type = 2      !  limiter type 
      muscl = 0         !  0 first order, 1 muscl ( high order )
      inmax = 8         !  maximum number of pressure correction steps
      maxout = 1000     !  maximum number of solver iteration
      omega = 1.0       !  relaxation parameter  ( 0.5 under - 2.0 over)
      implicit_fluxes   = .FALSE. ! fluxes are computed implicitly
      implicit_enthalpy = .FALSE. ! enthalpy solved implicitly

!
! :::::::::::::::::::::::  R E A D   N A M E L I S T S ::::::::::::::::
!

      IF( which == 'PP' ) THEN 
        IF( mpime == root ) THEN
          OPEN( UNIT=iuni_nml, FILE='pdac.nml', STATUS='UNKNOWN')
          OPEN( UNIT=iuni_grx, FILE='pdac.grx', STATUS='UNKNOWN')
          OPEN( UNIT=iuni_gry, FILE='pdac.gry', STATUS='UNKNOWN')
          OPEN( UNIT=iuni_grz, FILE='pdac.grz', STATUS='UNKNOWN')
          OPEN( UNIT=iuni_fld, FILE='pdac.fld', STATUS='UNKNOWN')
          WRITE(iuni_nml, * ) '<?xml version="1.0" encoding="UTF-8"?>'
          CALL iotk_write_begin( iuni_nml, "input" )
        END IF
      END IF
!
! ... Control Namelist ................................................
!
      IF(mpime == root) THEN

        READ(iunit, control) 

        IF( which == 'PP' ) THEN
          ! WRITE(iuni_nml, control) 
       
          CALL iotk_write_begin( iuni_nml, "control" )
            CALL iotk_write_begin( iuni_nml, "run_name" )
               WRITE( iuni_nml, * ) run_name
            CALL iotk_write_end( iuni_nml, "run_name" )
            CALL iotk_write_begin( iuni_nml, "job_type" )
               WRITE( iuni_nml, * ) job_type
            CALL iotk_write_end( iuni_nml, "job_type" )
            CALL iotk_write_begin( iuni_nml, "restart_mode" )
               WRITE( iuni_nml, * ) restart_mode
            CALL iotk_write_end( iuni_nml, "restart_mode" )
            CALL iotk_write_dat( iuni_nml, "time", time )
            CALL iotk_write_dat( iuni_nml, "tstop", tstop )
            CALL iotk_write_dat( iuni_nml, "dt", dt )
            CALL iotk_write_dat( iuni_nml, "lpr", lpr )
            CALL iotk_write_dat( iuni_nml, "tpr", tpr )
            CALL iotk_write_dat( iuni_nml, "tdump", tdump )
            CALL iotk_write_dat( iuni_nml, "nfil", irex )
            CALL iotk_write_dat( iuni_nml, "old_restart", old_restart )
            CALL iotk_write_dat( iuni_nml, "max_seconds", max_seconds )
          CALL iotk_write_end( iuni_nml, "control" )
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
      CALL bcast_logical(formatted_output,1,root)
      CALL bcast_logical(old_restart,1,root)
      CALL bcast_real(max_seconds,1,root)

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
        CASE DEFAULT
          CALL error(' input ',' unknown job_type '//TRIM(job_type), 1 )
      END SELECT

      CALL subsc_setup( job_type )
!
! ... Model Namelist ................................................
!
      IF(mpime == root) THEN

        READ(iunit, model) 

        IF( which == 'PP' ) THEN
          CALL iotk_write_begin( iuni_nml, "model" )
            CALL iotk_write_dat( iuni_nml, "gas_viscosity", gas_viscosity )
            CALL iotk_write_dat( iuni_nml, "part_viscosity", part_viscosity )
            CALL iotk_write_dat( iuni_nml, "iss", iss )
            CALL iotk_write_dat( iuni_nml, "repulsive_model",repulsive_model )
            CALL iotk_write_dat( iuni_nml, "iturb", iturb )
            CALL iotk_write_dat( iuni_nml, "modturbo", modturbo )
            CALL iotk_write_dat( iuni_nml, "cmut", cmut )
            CALL iotk_write_dat( iuni_nml, "rlim", rlim )
            CALL iotk_write_dat( iuni_nml, "gravx", gravx )
            CALL iotk_write_dat( iuni_nml, "gravz", gravz )
            CALL iotk_write_dat( iuni_nml, "ngas", ngas )
            CALL iotk_write_dat( iuni_nml, "density_specified", density_specified )
          CALL iotk_write_end( iuni_nml, "model" )
        END IF
      END IF

      CALL bcast_integer(irex,1,root)
      CALL bcast_integer(iss,1,root)
      CALL bcast_integer(repulsive_model,1,root)
      CALL bcast_integer(iturb,1,root)
      CALL bcast_integer(modturbo,1,root)
      CALL bcast_real(cmut,1,root)
      CALL bcast_real(rlim,1,root)
      CALL bcast_real(gravx,1,root)
      CALL bcast_real(gravy,1,root)
      CALL bcast_real(gravz,1,root)
      CALL bcast_integer(ngas,1,root)
      CALL bcast_logical(gas_viscosity,1,root)
      CALL bcast_logical(part_viscosity,1,root)
      CALL bcast_logical(density_specified,1,root)
!
! ... Mesh Namelist ...................................................
!
      IF(mpime == root) THEN
      
        READ(iunit, mesh) 
        
        IF( which == 'PP' ) THEN
          ! WRITE(iuni_nml, mesh) 
          CALL iotk_write_begin( iuni_nml, "mesh" )
            CALL iotk_write_dat( iuni_nml, "nx", nx )
            CALL iotk_write_dat( iuni_nml, "ny", ny )
            CALL iotk_write_dat( iuni_nml, "nz", nz )
            CALL iotk_write_dat( iuni_nml, "itc", itc )
            CALL iotk_write_dat( iuni_nml, "iuni", iuni )
            CALL iotk_write_dat( iuni_nml, "dx0", dx0 )
            CALL iotk_write_dat( iuni_nml, "dy0", dy0 )
            CALL iotk_write_dat( iuni_nml, "dz0", dz0 )
            CALL iotk_write_dat( iuni_nml, "origin_x", origin_x )
            CALL iotk_write_dat( iuni_nml, "origin_y", origin_y )
            CALL iotk_write_dat( iuni_nml, "origin_z", origin_z )
            CALL iotk_write_dat( iuni_nml, "mesh_partition", mesh_partition )
          CALL iotk_write_end( iuni_nml, "mesh" )
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
!
! ... Boundaries Namelist ................................................
!
      IF(mpime == root) THEN

        READ(iunit, boundaries) 

        IF( which == 'PP' ) THEN
          CALL iotk_write_begin( iuni_nml, "boundaries" )
          ! WRITE(iuni_nml, boundaries) 
            CALL iotk_write_dat( iuni_nml, "east", east )
            CALL iotk_write_dat( iuni_nml, "west", west )
            CALL iotk_write_dat( iuni_nml, "south", south )
            CALL iotk_write_dat( iuni_nml, "north", north )
            CALL iotk_write_dat( iuni_nml, "top", top )
            CALL iotk_write_dat( iuni_nml, "bottom", bottom )
            CALL iotk_write_dat( iuni_nml, "itp", itp )
            CALL iotk_write_begin( iuni_nml, "topography" )
              WRITE( iuni_nml, * ) topography
            CALL iotk_write_end( iuni_nml, "topography" )
            CALL iotk_write_dat( iuni_nml, "immb", immb )
          CALL iotk_write_end( iuni_nml, "boundaries" )
        END IF
      END IF

      CALL bcast_integer(east,1,root)
      CALL bcast_integer(west,1,root)
      CALL bcast_integer(north,1,root)
      CALL bcast_integer(south,1,root)
      CALL bcast_integer(bottom,1,root)
      CALL bcast_integer(top,1,root)
      CALL bcast_integer(itp,1,root)
      CALL bcast_character(topography,80,root)
      CALL bcast_integer(immb,1,root)
!
! ... Particles Namelist ..............................................
!
      IF(mpime == root) THEN

        READ(iunit, particles) 

        IF( which == 'PP' ) THEN
          ! WRITE(iuni_nml, particles) 
          CALL iotk_write_begin( iuni_nml, "particles" )
            CALL iotk_write_dat( iuni_nml, "nsolid", nsolid )
            CALL iotk_write_dat( iuni_nml, "diameter", diameter )
            CALL iotk_write_dat( iuni_nml, "density", density )
            CALL iotk_write_dat( iuni_nml, "sphericity", sphericity )
            CALL iotk_write_dat( iuni_nml, "viscosity", viscosity )
            CALL iotk_write_dat( iuni_nml, "specific_heat", specific_heat )
            CALL iotk_write_dat( iuni_nml, "thermal_conductivity", &
                                            thermal_conductivity )
          CALL iotk_write_end( iuni_nml, "particles" )
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
!
! ... Numeric Namelist ................................................
!
      IF(mpime == root) THEN

        READ(iunit, numeric) 

        IF( which == 'PP' ) THEN
          ! WRITE(iuni_nml, numeric) 
          CALL iotk_write_begin( iuni_nml, "numeric" )
            CALL iotk_write_dat( iuni_nml, "rungekut", rungekut )
            CALL iotk_write_dat( iuni_nml, "beta", beta )
            CALL iotk_write_dat( iuni_nml, "muscl", muscl )
            CALL iotk_write_dat( iuni_nml, "lim_type", lim_type )
            CALL iotk_write_dat( iuni_nml, "linmax", inmax )
            CALL iotk_write_dat( iuni_nml, "maxout", maxout )
            CALL iotk_write_dat( iuni_nml, "omega", omega )
            CALL iotk_write_dat( iuni_nml, "implicit_fluxes",   &
                                          & implicit_fluxes )
            CALL iotk_write_dat( iuni_nml, "implicit_enthalpy", &
                                          & implicit_enthalpy )
          CALL iotk_write_end( iuni_nml, "numeric" )
        END IF
      END IF

      CALL bcast_integer(rungekut,1,root)
      CALL bcast_real(beta,1,root)
      CALL bcast_integer(lim_type,1,root)
      CALL bcast_integer(muscl,1,root)
      CALL bcast_integer(inmax,1,root)
      CALL bcast_integer(maxout,1,root)
      CALL bcast_real(omega,1,root)
      CALL bcast_logical(implicit_fluxes,1,root)
      CALL bcast_logical(implicit_enthalpy,1,root)
!
! ... PP Namelist .....................................................
!
      IF( which == 'PP' ) THEN
        IF(mpime == root) THEN
          READ(iunit, pp) 
          ! WRITE(iuni_nml, pp) 
          CALL iotk_write_begin( iuni_nml, "pp" )
            CALL iotk_write_dat( iuni_nml, "first_out", first_out )
            CALL iotk_write_dat( iuni_nml, "last_out", last_out )
            CALL iotk_write_dat( iuni_nml, "incr_out", incr_out )
          CALL iotk_write_end( iuni_nml, "pp" )
        END IF
        CALL bcast_integer(first_out,1,root)
        CALL bcast_integer(last_out,1,root)
        CALL bcast_integer(incr_out,1,root)
      END IF
!
! :::::::::::::::::::::::  R E A D   C A R D S ::::::::::::::::::::::::
!
!
! ... Roughness Card ..................................................
!
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
        IF( which == 'PP' ) THEN
          CALL iotk_write_begin( iuni_nml, "roughness" )
            CALL iotk_write_dat( iuni_nml, "ir", zrough%ir )
            CALL iotk_write_dat( iuni_nml, "roucha", zrough%roucha )
            CALL iotk_write_dat( iuni_nml, "r", zrough%r( 1 : zrough%ir ) )
          CALL iotk_write_end( iuni_nml, "roughness" )
        END IF
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
! ... Mesh Card (cell sizes) ..........................................
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
!
! :::::::::::::::::::::   write post-processing files   :::::::::::::::::::::::
!
      IF(mpime == root) THEN

        IF( which == 'PP' ) THEN

          ALLOCATE( grx( nx - 1 ) )
          ALLOCATE( gry( ny - 1 ) )
          ALLOCATE( grz( nz - 1 ) )

          grx( 1 ) = origin_x
          WRITE( iuni_grx, * ) grx( 1 )
          DO i = 2, nx - 1
            grx( i ) = grx( i - 1 ) + delta_x( i )
            WRITE( iuni_grx, * ) grx( i )
          END DO
  
          gry( 1 ) = origin_y
          WRITE( iuni_gry, * ) gry( 1 )
          DO i = 2, ny - 1
            gry( i ) = gry( i - 1 ) + delta_y( i )
            WRITE( iuni_gry, * ) gry( i )
          END DO
  
          grz( 1 ) = origin_z
          WRITE( iuni_grz, * ) grz( 1 )
          DO i = 2, nz - 1
            grz( i ) = grz( i - 1 ) + delta_z( i )
            WRITE( iuni_grz, * ) grz( i )
          END DO

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

          CALL iotk_write_begin( iuni_nml, "mesh" )
            CALL iotk_write_dat( iuni_nml, "grx", grx )
            CALL iotk_write_dat( iuni_nml, "gry", gry )
            CALL iotk_write_dat( iuni_nml, "grz", grz )
          CALL iotk_write_end( iuni_nml, "mesh" )

          DEALLOCATE( grx )
          DEALLOCATE( gry )
          DEALLOCATE( grz )

        END IF

      END IF
!
! ... Fixed Flows Card (Inlet conditions) .............................
!
      npr = 0
      tend = .FALSE.

      IF(mpime == root) THEN

        fixed_flows_search: DO
          READ( 5, *, END = 300 ) card
          IF( TRIM(card) == 'FIXED_FLOWS' ) THEN
            EXIT fixed_flows_search
          END IF
        END DO fixed_flows_search

        READ(5,*) number_of_block

        IF (job_type == '2D') THEN

          DO n = 1, number_of_block
            READ(5,*) block_type(n), block_bounds(1,n),  block_bounds(2,n),  &
                      block_bounds(5,n), block_bounds(6,n)
            IF( block_type(n) == 1 .OR. block_type(n) == 5) THEN
              READ(5,*) fixed_vgas_x(n), fixed_vgas_z(n), fixed_pressure(n), &
                        fixed_gaseps(n), fixed_gastemp(n)
              READ(5,*) (fixed_vpart_x(k,n), fixed_vpart_z(k,n), &
                         fixed_parteps(k,n), fixed_parttemp(k,n), k=1, nsolid)
              READ(5,*) ( fixed_gasconc(ig,n), ig=1, max_ngas )
            ELSE IF( block_type(n) == 7) THEN
              CALL read_profile(n)
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
              READ(5,*) (fixed_vpart_x(k,n),fixed_vpart_y(k,n),fixed_vpart_z(k,n), &
                          fixed_parteps(k,n), fixed_parttemp(k,n), k=1, nsolid)
              READ(5,*) (fixed_gasconc(ig,n), ig=1, max_ngas )
            ENDIF
          END DO

        ELSE

          CALL error( ' input # FIXED FLOW ', ' unknown job_type', 1 )

        END IF

        IF( which == 'PP' ) THEN
          attr = ' '
          CALL iotk_write_attr( attr, "number_of_block", number_of_block )
          CALL iotk_write_begin( iuni_nml, "fixed_flows", attr )
          DO n = 1, number_of_block
            attr = ' '
            CALL iotk_write_attr( attr, "id", n )
            CALL iotk_write_attr( attr, "block_type", block_type(n) )
            CALL iotk_write_begin( iuni_nml, "block", attr )
              CALL iotk_write_dat( iuni_nml, "xlo", block_bounds(1,n) )
              CALL iotk_write_dat( iuni_nml, "xhi", block_bounds(2,n) )
              CALL iotk_write_dat( iuni_nml, "ylo", block_bounds(3,n) )
              CALL iotk_write_dat( iuni_nml, "yhi", block_bounds(4,n) )
              CALL iotk_write_dat( iuni_nml, "zlo", block_bounds(5,n) )
              CALL iotk_write_dat( iuni_nml, "zhi", block_bounds(6,n) )
              CALL iotk_write_dat( iuni_nml, "fixed_vgas_x", fixed_vgas_x(n) )
              CALL iotk_write_dat( iuni_nml, "fixed_vgas_y", fixed_vgas_y(n) )
              CALL iotk_write_dat( iuni_nml, "fixed_vgas_z", fixed_vgas_z(n) )
              CALL iotk_write_dat( iuni_nml, "fixed_pressure", fixed_pressure(n) )
              CALL iotk_write_dat( iuni_nml, "fixed_gaseps", fixed_gaseps(n) )
              CALL iotk_write_dat( iuni_nml, "fixed_gastemp", fixed_gastemp(n) )
              CALL iotk_write_dat( iuni_nml, "fixed_gasconc", fixed_gasconc( 1:max_ngas, n) )
              CALL iotk_write_dat( iuni_nml, "fixed_vpart_x", fixed_vpart_x( 1:nsolid, n) )
              CALL iotk_write_dat( iuni_nml, "fixed_vpart_y", fixed_vpart_y( 1:nsolid, n) )
              CALL iotk_write_dat( iuni_nml, "fixed_vpart_z", fixed_vpart_z( 1:nsolid, n) )
              CALL iotk_write_dat( iuni_nml, "fixed_parteps", fixed_parteps( 1:nsolid, n) )
              CALL iotk_write_dat( iuni_nml, "fixed_parttemp", fixed_parttemp( 1:nsolid, n) )
            CALL iotk_write_end( iuni_nml, "block" )
          END DO
          CALL iotk_write_end( iuni_nml, "fixed_flows" )
        END IF

        GOTO 310
 300    tend = .TRUE.
 310    CONTINUE

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
 
      CALL bcast_integer(npr, 1, root)
      CALL bcast_real(prof_vgas_x, SIZE(prof_vgas_x), root)
      CALL bcast_real(prof_vgas_z, SIZE(prof_vgas_z), root)
      CALL bcast_real(prof_pressure, SIZE(prof_pressure), root)
      CALL bcast_real(prof_gaseps, SIZE(prof_gaseps), root)
      CALL bcast_real(prof_gastemp, SIZE(prof_gastemp), root)
      CALL bcast_real(prof_vpart_x, SIZE(prof_vpart_x), root)
      CALL bcast_real(prof_vpart_z, SIZE(prof_vpart_z), root)
      CALL bcast_real(prof_parteps, SIZE(prof_parteps), root)
      CALL bcast_real(prof_parttemp, SIZE(prof_parttemp), root)
      CALL bcast_real(prof_gasconc, SIZE(prof_gasconc), root)
!
! ... Initial conditions Card .........................................
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
          READ(5,*) initial_vgas_x, initial_vgas_z, initial_pressure, &
            initial_void_fraction, max_packing, initial_temperature
          READ(5,*) initial_vpart_x, initial_vpart_z 
          READ(5,*) (initial_gasconc(ig), ig=1, max_ngas)
        ELSE IF( job_type == '3D' ) THEN
          READ(5,*) initial_vgas_x, initial_vgas_y, initial_vgas_z, &
            initial_pressure, initial_void_fraction, max_packing, initial_temperature
          READ(5,*) initial_vpart_x, initial_vpart_y, initial_vpart_z
          READ(5,*) (initial_gasconc(ig), ig=1, max_ngas)
        ELSE 
          CALL error('input # INITIAL_CONDITIONS', 'unknown job_type',1) 
        ENDIF

        IF( which == 'PP' ) THEN
          CALL iotk_write_begin( iuni_nml, "initial_conditions" )
            CALL iotk_write_dat( iuni_nml, "initial_vgas_x", initial_vgas_x )
            CALL iotk_write_dat( iuni_nml, "initial_vgas_y", initial_vgas_y )
            CALL iotk_write_dat( iuni_nml, "initial_vgas_z", initial_vgas_z )
            CALL iotk_write_dat( iuni_nml, "initial_pressure", initial_pressure )
            CALL iotk_write_dat( iuni_nml, "initial_void_fraction", initial_void_fraction )
            CALL iotk_write_dat( iuni_nml, "max_packing", max_packing )
            CALL iotk_write_dat( iuni_nml, "initial_temperature", initial_temperature )
            CALL iotk_write_dat( iuni_nml, "initial_vpart_x", initial_vpart_x )
            CALL iotk_write_dat( iuni_nml, "initial_vpart_y", initial_vpart_y )
            CALL iotk_write_dat( iuni_nml, "initial_vpart_z", initial_vpart_z )
            CALL iotk_write_dat( iuni_nml, "initial_gasconc", initial_gasconc( 1 : max_ngas ) )
          CALL iotk_write_end( iuni_nml, "initial_conditions" )
        END IF

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
          CALL iotk_write_end( iuni_nml, "input" )
          CLOSE( UNIT=iuni_nml )
          CLOSE( UNIT=iuni_fld )
          CLOSE( UNIT=iuni_grx )
          CLOSE( UNIT=iuni_gry )
          CLOSE( UNIT=iuni_grz )
        END IF
      END IF

      END SUBROUTINE input
!----------------------------------------------------------------------
      SUBROUTINE initc

      USE atmosphere, ONLY: v0, u0, w0, p0, temp0, us0, vs0, ws0, &
     &                      ep0, epsmx0
      USE control_flags, ONLY: job_type, lpr
      USE dimensions
      USE grid, ONLY: dx, dy, dz, itc
      USE grid, ONLY: iob, zzero
      USE initial_conditions, ONLY: epsob, tpob, ygc0, ygcob,   &
     &     ugob, vgob, wgob, upob, vpob, wpob, pob, tgob, epob
      USE initial_conditions, ONLY: ugpr, wgpr, ppr, eppr, tgpr, &
                                    uppr,wppr,epspr,tppr, ygcpr, npr
      USE particles_constants, ONLY: rl, inrl, kap, cmus, phis, cps, dk

      IMPLICIT NONE

      INTEGER :: ig, is

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
      ygcob(1:max_ngas,1:no) = fixed_gasconc(1:max_ngas,1:no)

      ugpr(1:npr)  = prof_vgas_x(1:npr)
      wgpr(1:npr)  = prof_vgas_z(1:npr)
      ppr(1:npr)  = prof_pressure(1:npr)
      eppr(1:npr)  = prof_gaseps(1:npr)
      tgpr(1:npr)  = prof_gastemp(1:npr)
      uppr(1:nsolid,1:npr) = prof_vpart_x(1:nsolid,1:npr)
      wppr(1:nsolid,1:npr) = prof_vpart_z(1:nsolid,1:npr)
      epspr(1:nsolid,1:npr) = prof_parteps(1:nsolid,1:npr)
      tppr(1:nsolid,1:npr) = prof_parttemp(1:nsolid,1:npr)
      ygcpr(1:max_ngas,1:npr) = prof_gasconc(1:max_ngas,1:npr)

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

      ygc0(1:max_ngas) = initial_gasconc(1:max_ngas)
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
      SUBROUTINE read_profile(n)

      USE dimensions, ONLY: max_ngas, nsolid
      USE initial_conditions, ONLY: npr
      USE parallel, ONLY: mpime, root

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      INTEGER :: k, ig, is
      INTEGER :: x1,x2,z1,z2

      OPEN(UNIT=17, FILE='input.prof', STATUS='OLD')
      READ(17,*) x1,x2,z1,z2

      npr = z2 - z1 + 1

      IF ( (z2-z1) /= (block_bounds(6,n)-block_bounds(5,n)) .OR. (x1/=x2) ) THEN
        CALL error('setup','Error in input profile, block:', n)
      END IF

      DO k= 1, npr
        READ(17,*) prof_vgas_x(k),prof_vgas_z(k),prof_pressure(k),  &
                   prof_gaseps(k),prof_gastemp(k),                  &
                  (prof_vpart_x(is,k),prof_vpart_z(is,k),           &
                   prof_parteps(is,k),prof_parttemp(is,k), is=1,nsolid), &
                  (prof_gasconc(ig,k), ig = 1,max_ngas)
      END DO
      CLOSE(17)

      END SUBROUTINE read_profile
!----------------------------------------------------------------------
      END MODULE input_module
!----------------------------------------------------------------------
