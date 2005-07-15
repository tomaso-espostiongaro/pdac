!----------------------------------------------------------------------
      MODULE input_module
!----------------------------------------------------------------------
      USE dimensions, ONLY: max_nsolid, ngas, nroughx, max_size, &
          max_nblock, max_ngas, nz, nx, ny
      USE io_files, ONLY: logunit

      PRIVATE
      PUBLIC :: input, initc, number_of_block, run_name
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
      INTEGER :: mass_limiter, vel_limiter

! ... MESH
      REAL*8 :: delta_x(max_size)
      REAL*8 :: delta_y(max_size)
      REAL*8 :: delta_z(max_size)

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
!
      REAL*8 :: atm_O2
      REAL*8 :: atm_N2
      REAL*8 :: atm_CO2
      REAL*8 :: atm_H2
      REAL*8 :: atm_H2O
      REAL*8 :: atm_Air
      REAL*8 :: atm_SO2
!
      REAL*8 :: troposphere_z
      REAL*8 :: tropopause_z
      REAL*8 :: lower_stratosphere_z
      REAL*8 :: upper_stratosphere_z
      REAL*8 :: ozone_layer_z
      REAL*8 :: lower_mesosphere_z
      REAL*8 :: upper_mesosphere_z
!
      REAL*8 :: troposphere_grad
      REAL*8 :: tropopause_grad
      REAL*8 :: lower_stratosphere_grad
      REAL*8 :: upper_stratosphere_grad
      REAL*8 :: ozone_layer_grad
      REAL*8 :: lower_mesosphere_grad
      REAL*8 :: upper_mesosphere_grad
!
      REAL*8 :: vent_O2
      REAL*8 :: vent_N2
      REAL*8 :: vent_CO2
      REAL*8 :: vent_H2
      REAL*8 :: vent_H2O
      REAL*8 :: vent_Air
      REAL*8 :: vent_SO2
!
      REAL*8 :: dome_O2
      REAL*8 :: dome_N2
      REAL*8 :: dome_CO2
      REAL*8 :: dome_H2
      REAL*8 :: dome_H2O
      REAL*8 :: dome_Air
      REAL*8 :: dome_SO2
!
      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE input( iunit )

      USE atmospheric_conditions, ONLY: gravx, gravy, gravz, stratification
      USE atmospheric_conditions, ONLY: wind_x, wind_y, wind_z, &
          p_ground, t_ground, void_fraction, max_packing
      USE blunt_body, ONLY: ibl, nblu
      USE control_flags, ONLY: job_type, lpr, imr
      USE control_flags, ONLY: implicit_fluxes, implicit_enthalpy
      USE domain_decomposition, ONLY: mesh_partition
      USE dome_conditions, ONLY: xdome, ydome, dome_volume, temperature, particle_fraction, overpressure, &
          idome, gas_flux, permeability, dome_gasvisc, idw
      USE enthalpy_matrix, ONLY: flim, tforce
      USE eos_gas, ONLY: update_eosg
      USE flux_limiters, ONLY: beta, muscl
      USE gas_solid_viscosity, ONLY: gas_viscosity, part_viscosity
      USE gas_solid_viscosity, ONLY: repulsive_model
      USE grid, ONLY: dx, dy, dz, itc, zzero
      USE grid, ONLY: west, east, south, north, bottom, top
      USE grid, ONLY: domain_x, domain_y, domain_z, n0x, n0y, n0z
      USE grid, ONLY: alpha_x, alpha_y, alpha_z
      USE grid, ONLY: center_x, center_y
      USE grid, ONLY: dxmin, dxmax, dymin, dymax, dzmin, dzmax
      USE grid, ONLY: maxbeta, grigen
      USE initial_conditions, ONLY: density_specified
      USE vent_conditions, ONLY: u_gas,v_gas,w_gas,p_gas,t_gas, wrat, &
          u_solid, v_solid, w_solid,  ep_solid, t_solid, base_radius, &
          crater_radius, vent_radius, xvent, yvent, ivent, iali, irand, &
          ipro, rad_file
      USE immersed_boundaries, ONLY: immb
      USE iterative_solver, ONLY: inmax, maxout, omega, optimization, delg, vforce
      USE io_restart, ONLY: max_seconds, nfil
      USE output_dump, ONLY: formatted_output
      USE parallel, ONLY: mpime, root
      USE particles_constants, ONLY: rl, inrl, kap, &
     &     cmus, phis, cps, dk, nsolid
      USE phases_matrix, ONLY: rlim
      USE reactions, ONLY: irex
      USE roughness_module, ONLY: zrough
      USE set_indexes, ONLY: subsc_setup
      USE specific_heat_module, ONLY: icpc
      USE time_parameters, ONLY: time, tstop, dt, tpr, tdump, itd, & 
     &                            timestart, rungekut, tau
      USE turbulence_model, ONLY: iturb, cmut, iss, modturbo
      USE volcano_topography, ONLY: itp, iavv, cellsize, filtersize, &
                         dem_file, nocrater, rim_quota, ismt, itrans
!
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: iunit

      NAMELIST / control / run_name, job_type, restart_mode,       &
        time, tstop, dt, lpr, imr, tpr, tdump, nfil, tau,          &
        formatted_output, max_seconds

      NAMELIST / model / icpc, irex, gas_viscosity, part_viscosity,      &
        iss, repulsive_model, iturb, modturbo, cmut,                     &
        gravx, gravy, gravz, ngas, density_specified

      NAMELIST / mesh / nx, ny, nz, itc, iuni, dx0, dy0, dz0, zzero, &
        center_x, center_y, alpha_x, alpha_y, alpha_z,  &
        dxmin, dxmax, dymin, dymax, dzmin, dzmax, n0x, n0y, n0z,     &
        domain_x, domain_y, domain_z, maxbeta, grigen, mesh_partition

      NAMELIST / boundaries / west, east, south, north, bottom, top, &
        immb, ibl

      NAMELIST / topography / dem_file, itp, iavv, nocrater, &
        rim_quota, filtersize, cellsize, ismt, zrough, itrans
      
      NAMELIST / inlet / ivent, iali, irand, ipro, rad_file, wrat, &
        crater_radius, &
        xvent, yvent, vent_radius, base_radius, u_gas, v_gas, w_gas,  &
        p_gas, t_gas, u_solid, v_solid, w_solid, ep_solid, t_solid, &
        vent_O2, vent_N2, vent_CO2, vent_H2, vent_H2O, vent_Air, vent_SO2

      NAMELIST / dome / xdome, ydome, dome_volume, temperature, particle_fraction, idome, overpressure, &
                        gas_flux, permeability, dome_gasvisc, idw

      NAMELIST / atmosphere / wind_x, wind_y, wind_z, p_ground, t_ground, &
        void_fraction, max_packing, atm_O2, atm_N2, atm_CO2, atm_H2, atm_H2O, &
        atm_Air, atm_SO2, troposphere_grad, tropopause_grad, &
        lower_stratosphere_grad, upper_stratosphere_grad, ozone_layer_grad, &
        lower_mesosphere_grad, upper_mesosphere_grad, troposphere_z, &
        tropopause_z, lower_stratosphere_z, upper_stratosphere_z, &
        ozone_layer_z, lower_mesosphere_z, upper_mesosphere_z, stratification

      NAMELIST / particles / nsolid, diameter, density, sphericity, &
        viscosity, specific_heat, thermal_conductivity

      NAMELIST / numeric / rungekut, beta, muscl, mass_limiter, vel_limiter, &
        inmax, maxout, omega, delg, implicit_fluxes, implicit_enthalpy, &
        update_eosg, optimization, lim_type, tforce, rlim, flim, vforce

      INTEGER :: i, j, k, n, m, ig, ierr, lim_type
      CHARACTER(LEN=80) :: card
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
      lpr = 0           ! verbosity 
      imr = 0           ! =1 print mass residuals
      tpr = 1.0D0       ! write to output file every tpr seconds of simulated time
      tdump = 20.0D0    ! write restart every tdump seconds of simulated time
      nfil = 0          ! output file index
      formatted_output = .TRUE.
      max_seconds = 40000.0
      tau = 0.D0

! ... Model

      icpc = 1          ! ( 1 specific heat depends on temperature)
      irex = 1          ! ( 1 no reaction, 2 use hrex. Not used )
      gas_viscosity  = .TRUE. ! include molecular gas viscosity
      part_viscosity = .TRUE. ! include collisional particle viscosity
      density_specified = .FALSE. ! density specified instead of temperature
      iss  = 0          ! ( 0 no solid turbulent stress, 1 sub-grid stress) 
      repulsive_model = 1 ! ( 0 no Coulombic repulsive model )
      iturb = 1         ! turbulence  ( 0 no turbo, 1 turbo, 2 turbo + rough )
      modturbo = 1      ! turbulence  ( 1 smag, 2 dynamic )
      cmut = 0.1D0      ! Smagorinsky constant
      gravx = 0.0D0     ! gravity along x
      gravy = 0.0D0     ! gravity along y
      gravz = -9.81D0   ! gravity along z
      ngas = 2          ! max number of gas components

! ... Mesh

      nx = 100                !  number of cell in the X directions
      nz = 100                !  number of cell in the Z directions
      ny = 1                  !  number of cell in the Y directions
      n0x = 10                !  number of cell with minimum size
      n0y = 10                !  number of cell with minimum size
      n0z = 10                !  number of cell with minimum size
      itc = 0                 !  itc = 1 cylindrical coordinates are used
      grigen = 0              !  flag for grid generation 0 = no grid gen.
      maxbeta = 1.2           !  maximum increase rate for non-uniform meshes
      mesh_partition = 1      !  type of partition
      iuni = 0                !  1 = uniform grid, 0 = non uniform grid
      dz0  = 20.D0            !  default cell z size in meters
      dx0  = 20.D0            !  default cell x size in meters
      dy0  = 20.D0            !  default cell y size in meters
      dxmin = 10.D0           !  default cell x minimum size in metres
      dymin = 10.D0           !  default cell y minimum size in metres
      dzmin = 10.D0           !  default cell z minimum size in metres
      dxmax = 100.D0           !  default cell x maximum size in metres
      dymax = 100.D0           !  default cell y maximum size in metres
      dzmax = 100.D0           !  default cell z maximum size in metres
      alpha_x  = 0.5D0        !  default x coo. of the mesh center
      alpha_y  = 0.5D0        !  default y coo. of the mesh center
      alpha_z  = 0.0D0        !  default z coo. of the mesh center
      center_x  = 0.D0        !  default x coo. of the mesh center
      center_y  = 0.D0        !  default y coo. of the mesh center
      domain_x  = 10.D3       !  domain size in metres
      domain_y  = 10.D3       !  domain size in metres
      domain_z  = 10.D3       !  domain size in metres
      zzero = 0.0D0           !  grid bottom level

! ... Boundaries

      west = 2                !
      east = 2                !
      bottom = 3              ! boundary types
      top = 4                 !
      south = 2               !
      north = 2               !
      immb  = 0               ! 1: use immersed boundaries
      ibl  = 0                ! 1: compute drag and lift on blocks

! ... Topography

      dem_file = 'topo.dat' ! file containing the topographic profile
      itp   = 0               ! itp = 1 => read topography from file
      ismt   = 1               ! 1 = median filter; 2 = gaussian filter
      iavv   = 0              ! iavv = 1 => average volcano topography
      nocrater = .FALSE. ! flatten the crater
      itrans = .TRUE.         ! translate vertically
      rim_quota = 1000.D0     ! index of the rim quota 
      filtersize  = 50        ! low-pass filter size
      cellsize  = 10          ! resolution of the resized dem
      zrough = 1.D0           ! average roughness length

! ... Inlet

      ivent = 0               ! 0: specify inlet blocks 1: circular vent
      iali  = 0               ! 1: vent antialiasing ON
      irand = 0               ! 1: circular vent specified on average
      ipro = 0                ! 1: inlet radial profile
      rad_file = 'profile.rad'! file with the radial profile
      xvent  = 0.D0           ! coordinates of the vent
      yvent  = 0.D0           ! coordinates of the vent
      vent_radius = 100.D0    ! vent radius
      base_radius = 200.D0    ! base radius
      crater_radius = 500.D0  ! maximum radius of the external crater rim
      u_gas = 0.D0            ! gas velocity x
      v_gas = 0.D0            ! gas velocity y
      w_gas = 0.D0            ! gas velocity z
      wrat = 1.D0             ! maximum/mean vertical velocity ratio
      p_gas = 1.01325D5       ! gas pressure
      t_gas  = 288.15D0       ! gas temperature
      u_solid = 0.D0           ! particle velocity x (array)
      v_solid = 0.D0           ! particle velocity y (array)
      w_solid = 0.D0           ! particle velocity z (array)
      ep_solid = 0.D0          ! particle fraction   (array)
      t_solid  = 288.15D0      ! gas temperature
      vent_O2  = 0.D0          ! gas components mass fractions
      vent_N2  = 0.D0
      vent_CO2 = 0.D0
      vent_H2  = 0.D0
      vent_H2O = 1.D0
      vent_Air = 0.D0
      vent_SO2 = 0.D0

! ... Dome

      idome = 0               ! Flag for automatic dome conditions
      idw = 0                 ! Flag for adding dome hydrostatic pressure
      xdome = 0.0             ! UTM longitude of the dome center
      ydome = 0.0             ! UTM latitude of the dome center
      dome_volume = 0.0        ! total volume of exploded mass
      overpressure = 100.D5             ! overpressure of the dome
      particle_fraction = 0.D0          ! particle fractions
      gas_flux = 400.D0        ! gas flux through the conduit
      temperature = 1100.D0        ! gas temperature
      permeability = 1.D-12    ! permeability of the dome
      dome_gasvisc = 1.D-5        ! viscosity of the gas
      dome_O2  = 0.D0          ! gas components mass fractions
      dome_N2  = 0.D0
      dome_CO2 = 0.D0
      dome_H2  = 0.D0
      dome_H2O = 1.D0
      dome_Air = 0.D0
      dome_SO2 = 0.D0

! ... Atmosphere

      wind_x = 0.D0            ! wind velocity x
      wind_y = 0.D0            ! wind velocity y
      wind_z = 0.D0            ! wind velocity z
      stratification = .TRUE.  ! set atmospheric stratification
      p_ground  = 1.01325D5    ! atmospheric pressure at ground level
      t_ground  = 288.15D0     ! atmospheric temperature at ground level
      void_fraction = 1.D0     ! void fraction
      max_packing = 0.6413     ! maximum packing vol. fract.
      atm_O2  = 0.D0           ! gas components mass fractions
      atm_N2  = 0.D0
      atm_CO2 = 0.D0
      atm_H2  = 0.D0
      atm_H2O = 0.D0
      atm_Air = 1.D0
      atm_SO2 = 0.D0
      troposphere_z = 1.1D4    ! top level of atmospheric layers
      tropopause_z = 2.0D4 
      lower_stratosphere_z = 3.2D4
      upper_stratosphere_z = 4.7D4
      ozone_layer_z = 5.1D4
      lower_mesosphere_z = 7.1D4
      upper_mesosphere_z = 8.0D4
      troposphere_grad = -6.5D-3  ! temperature gradients in atmosph. layers
      tropopause_grad = 0.D0
      lower_stratosphere_grad = 1.D-3
      upper_stratosphere_grad = 2.8D-3
      ozone_layer_grad = 0.D0
      lower_mesosphere_grad = -2.8D-3 
      upper_mesosphere_grad = -2.0D-3 

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
      mass_limiter = 2  !  limiter type in mass equation
      vel_limiter = 2   !  limiter type in momentum equation
      muscl = 0         !  0 first order, 1 muscl ( high order )
      inmax = 8         !  maximum number of pressure correction steps
      maxout = 500      !  maximum number of solver iteration
      delg = 1.D-8      !  residual limit relative to gas bulk density
      omega = 1.1       !  relaxation parameter  ( 0.5 under - 2.0 over)
      optimization = 1  !  optimization degree on iterative solver
      tforce  = .FALSE. !  force temperatures
      vforce  = .FALSE. !  force velocities
      implicit_fluxes   = .FALSE. ! fluxes are computed implicitly
      implicit_enthalpy = .FALSE. ! enthalpy solved implicitly
      update_eosg       = .FALSE.  ! update density after temperature
      rlim = 1.0D-8     ! 
                        ! limit for off-diagonal contribution in matrix
      flim = 1.0D-6     ! inversion
                        ! 

!

      IF(mpime == root) WRITE(logunit,*) 'Reading input ...'
!
! :::::::::::::::::::::::  R E A D   N A M E L I S T S ::::::::::::::::
!
!
! ... Control Namelist ................................................
!
      IF(mpime == root) READ(iunit, control) 
!
      CALL bcast_character(run_name,80,root)
      CALL bcast_character(job_type,80,root)
      CALL bcast_character(restart_mode,80,root)
      CALL bcast_real(time,1,root)
      CALL bcast_real(tstop,1,root)
      CALL bcast_real(dt,1,root)
      CALL bcast_integer(lpr,1,root)
      CALL bcast_integer(imr,1,root)
      CALL bcast_real(tpr,1,root)
      CALL bcast_real(tdump,1,root)
      CALL bcast_integer(nfil,1,root)
      CALL bcast_logical(formatted_output,1,root)
      CALL bcast_real(max_seconds,1,root)
      CALL bcast_real(tau,1,root)

      SELECT CASE ( TRIM(restart_mode) )
        CASE ('check_geom')
          itd = -1
        CASE ('check_init')
          itd = 0
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
      IF(mpime == root) READ(iunit, model) 

      CALL bcast_integer(icpc,1,root)
      CALL bcast_integer(irex,1,root)
      CALL bcast_integer(iss,1,root)
      CALL bcast_integer(repulsive_model,1,root)
      CALL bcast_integer(iturb,1,root)
      CALL bcast_integer(modturbo,1,root)
      CALL bcast_real(cmut,1,root)
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
      IF(mpime == root) READ(iunit, mesh) 

      CALL bcast_integer(nx,1,root)
      CALL bcast_integer(ny,1,root)
      CALL bcast_integer(nz,1,root)
      CALL bcast_integer(n0x,1,root)
      CALL bcast_integer(n0y,1,root)
      CALL bcast_integer(n0z,1,root)
      CALL bcast_integer(itc,1,root)
      CALL bcast_integer(iuni,1,root)
      CALL bcast_integer(grigen,1,root)
      CALL bcast_integer(mesh_partition,1,root)
      CALL bcast_real(dx0,1,root)
      CALL bcast_real(dy0,1,root)
      CALL bcast_real(dz0,1,root)
      CALL bcast_real(dxmin,1,root)
      CALL bcast_real(dymin,1,root)
      CALL bcast_real(dzmin,1,root)
      CALL bcast_real(dxmax,1,root)
      CALL bcast_real(dymax,1,root)
      CALL bcast_real(dzmax,1,root)
      CALL bcast_real(alpha_x,1,root)
      CALL bcast_real(alpha_y,1,root)
      CALL bcast_real(alpha_z,1,root)
      CALL bcast_real(center_x,1,root)
      CALL bcast_real(center_y,1,root)
      CALL bcast_real(domain_x,1,root)
      CALL bcast_real(domain_y,1,root)
      CALL bcast_real(domain_z,1,root)
      CALL bcast_real(maxbeta,1,root)
      CALL bcast_real(zzero,1,root)
!
! ... Boundaries Namelist ................................................
!
      IF(mpime == root) READ(iunit, boundaries) 

      CALL bcast_integer(east,1,root)
      CALL bcast_integer(west,1,root)
      CALL bcast_integer(north,1,root)
      CALL bcast_integer(south,1,root)
      CALL bcast_integer(bottom,1,root)
      CALL bcast_integer(top,1,root)
      CALL bcast_integer(immb,1,root)
      CALL bcast_integer(ibl,1,root)
!
! ... Topography Namelist ................................................
!
      IF(mpime == root) READ(iunit, topography) 

      CALL bcast_character(dem_file,80,root)
      CALL bcast_integer(itp,1,root)
      CALL bcast_integer(ismt,1,root)
      CALL bcast_integer(iavv,1,root)
      CALL bcast_logical(nocrater,1,root)
      CALL bcast_logical(itrans,1,root)
      CALL bcast_real(rim_quota,1,root)
      CALL bcast_real(filtersize,1,root)
      CALL bcast_real(cellsize,1,root)
      CALL bcast_real(zrough,1,root)
!
! ... Inlet Namelist ................................................
!
      IF(mpime == root) READ(iunit, inlet) 

      CALL bcast_integer(ivent,1,root)
      CALL bcast_integer(iali,1,root)
      CALL bcast_integer(irand,1,root)
      CALL bcast_integer(ipro,1,root)
      CALL bcast_character(rad_file,80,root)
      CALL bcast_real(xvent,1,root)
      CALL bcast_real(yvent,1,root)
      CALL bcast_real(vent_radius,1,root)
      CALL bcast_real(base_radius,1,root)
      CALL bcast_real(crater_radius,1,root)
      CALL bcast_real(wrat,1,root)
      CALL bcast_real(u_gas,1,root)
      CALL bcast_real(v_gas,1,root)
      CALL bcast_real(w_gas,1,root)
      CALL bcast_real(p_gas,1,root)
      CALL bcast_real(t_gas,1,root)
      CALL bcast_real(u_solid,max_nsolid,root)
      CALL bcast_real(v_solid,max_nsolid,root)
      CALL bcast_real(w_solid,max_nsolid,root)
      CALL bcast_real(ep_solid,max_nsolid,root)
      CALL bcast_real(t_solid,max_nsolid,root)
      CALL bcast_real(vent_O2,1,root)
      CALL bcast_real(vent_N2,1,root)
      CALL bcast_real(vent_CO2,1,root)
      CALL bcast_real(vent_H2,1,root)
      CALL bcast_real(vent_H2O,1,root)
      CALL bcast_real(vent_Air,1,root)
      CALL bcast_real(vent_SO2,1,root)
!
! ... Inlet Namelist ................................................
!
      IF(mpime == root) READ(iunit, dome) 

      CALL bcast_integer(idome,1,root)
      CALL bcast_integer(idw,1,root)
      CALL bcast_real(xdome,1,root)
      CALL bcast_real(ydome,1,root)
      CALL bcast_real(dome_volume,1,root)
      CALL bcast_real(temperature,1,root)
      CALL bcast_real(overpressure,1,root)
      CALL bcast_real(particle_fraction,max_nsolid,root)
      CALL bcast_real(gas_flux,1,root)
      CALL bcast_real(permeability,1,root)
      CALL bcast_real(dome_gasvisc,1,root)
      CALL bcast_real(dome_O2,1,root)
      CALL bcast_real(dome_N2,1,root)
      CALL bcast_real(dome_CO2,1,root)
      CALL bcast_real(dome_H2,1,root)
      CALL bcast_real(dome_H2O,1,root)
      CALL bcast_real(dome_Air,1,root)
      CALL bcast_real(dome_SO2,1,root)
!
! ... Atmosphere Namelist ................................................
!
      IF(mpime == root) READ(iunit, atmosphere) 

      CALL bcast_real(wind_x,1,root)
      CALL bcast_real(wind_y,1,root)
      CALL bcast_real(wind_z,1,root)
      CALL bcast_logical(stratification,1,root)
      CALL bcast_real(p_ground,1,root)
      CALL bcast_real(t_ground,1,root)
      CALL bcast_real(void_fraction,1,root)
      CALL bcast_real(max_packing,1,root)
      CALL bcast_real(atm_O2,1,root)
      CALL bcast_real(atm_N2,1,root)
      CALL bcast_real(atm_CO2,1,root)
      CALL bcast_real(atm_H2,1,root)
      CALL bcast_real(atm_H2O,1,root)
      CALL bcast_real(atm_Air,1,root)
      CALL bcast_real(atm_SO2,1,root)
      CALL bcast_real(troposphere_z,1,root)
      CALL bcast_real(tropopause_z,1,root)
      CALL bcast_real(lower_stratosphere_z,1,root)
      CALL bcast_real(upper_stratosphere_z,1,root)
      CALL bcast_real(ozone_layer_z,1,root)
      CALL bcast_real(lower_mesosphere_z,1,root)
      CALL bcast_real(upper_mesosphere_z,1,root)
      CALL bcast_real(troposphere_grad,1,root)
      CALL bcast_real(tropopause_grad,1,root)
      CALL bcast_real(lower_stratosphere_grad,1,root)
      CALL bcast_real(upper_stratosphere_grad,1,root)
      CALL bcast_real(ozone_layer_grad,1,root)
      CALL bcast_real(lower_mesosphere_grad,1,root)
      CALL bcast_real(upper_mesosphere_grad,1,root)
!
! ... Particles Namelist ..............................................
!
      IF(mpime == root) READ(iunit, particles) 

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
      IF(mpime == root) READ(iunit, numeric) 

      CALL bcast_integer(rungekut,1,root)
      CALL bcast_real(beta,1,root)
      CALL bcast_integer(mass_limiter,1,root)
      CALL bcast_integer(vel_limiter,1,root)
      CALL bcast_integer(muscl,1,root)
      CALL bcast_integer(inmax,1,root)
      CALL bcast_integer(delg,1,root)
      CALL bcast_integer(maxout,1,root)
      CALL bcast_integer(optimization,1,root)
      CALL bcast_real(omega,1,root)
      CALL bcast_logical(tforce,1,root)
      CALL bcast_logical(vforce,1,root)
      CALL bcast_logical(implicit_fluxes,1,root)
      CALL bcast_logical(implicit_enthalpy,1,root)
      CALL bcast_logical(update_eosg,1,root)
      CALL bcast_real(rlim,1,root)
      CALL bcast_real(flim,1,root)

!
! :::::::::::::::::::::::  R E A D   C A R D S ::::::::::::::::::::::::
!
! ... Mesh Card (cell sizes) ..........................................
!
      tend = .FALSE.
      IF(mpime == root) THEN
        mesh_search: DO
          READ(iunit,*,END=200) card
          IF( TRIM(card) == 'MESH' ) THEN
            EXIT mesh_search
          END IF
        END DO mesh_search

        IF (iuni == 0) THEN
          IF (grigen == 0) THEN
            READ(iunit,*) (delta_x(i),i=1,nx)
            IF( job_type == '3D' ) THEN
              READ(iunit,*) (delta_y(j),j=1,ny)
            END IF
            READ(iunit,*) (delta_z(k),k=1,nz)
          END IF
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
! ... Fixed Flows Card (Inlet conditions) .............................
!
      tend = .FALSE.

      IF(mpime == root) THEN

        fixed_flows_search: DO
          READ( iunit, *, END = 300 ) card
          IF( TRIM(card) == 'FIXED_FLOWS' ) THEN
            EXIT fixed_flows_search
          END IF
        END DO fixed_flows_search

        READ(iunit,*) number_of_block
        IF (ibl >= 1) READ(iunit,*) nblu(1:number_of_block)

        IF (job_type == '2D') THEN

          DO n = 1, number_of_block
            READ(iunit,*) block_type(n), block_bounds(1,n),  block_bounds(2,n),  &
                      block_bounds(5,n), block_bounds(6,n)
            IF( block_type(n) == 1 .OR. block_type(n) == 5) THEN
              READ(iunit,*) fixed_vgas_x(n), fixed_vgas_z(n), fixed_pressure(n), &
                        fixed_gaseps(n), fixed_gastemp(n)
              READ(iunit,*) (fixed_vpart_x(k,n), fixed_vpart_z(k,n), &
                         fixed_parteps(k,n), fixed_parttemp(k,n), k=1, nsolid)
              READ(iunit,*) ( fixed_gasconc(ig,n), ig=1, max_ngas )
            ENDIF
          END DO

        ELSE IF (job_type == '3D') THEN

          DO n = 1, number_of_block
            READ(iunit,*) block_type(n), block_bounds(1,n), block_bounds(2,n), &
                                     block_bounds(3,n), block_bounds(4,n), &
                                     block_bounds(5,n), block_bounds(6,n)
            IF( block_type(n) == 1 .OR. block_type(n) == 5) THEN
              READ(iunit,*) fixed_vgas_x(n), fixed_vgas_y(n), fixed_vgas_z(n), &
                        fixed_pressure(n), fixed_gaseps(n), fixed_gastemp(n)
              READ(iunit,*) (fixed_vpart_x(k,n),fixed_vpart_y(k,n),fixed_vpart_z(k,n), &
                          fixed_parteps(k,n), fixed_parttemp(k,n), k=1, nsolid)
              READ(iunit,*) (fixed_gasconc(ig,n), ig=1, max_ngas )
            ENDIF
          END DO

        ELSE

          CALL error( ' input # FIXED FLOW ', ' unknown job_type', 1 )

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
      IF (number_of_block > 0) THEN
        CALL bcast_integer(nblu, SIZE(nblu), root)
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
      END IF
!.....................................................................
!
      IF (lpr > 1) CALL write_xml


      IF(mpime == root) &
        WRITE(logunit,*) 'END input'


!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE write_xml
      !
      ! ... Writes input parameters into an XML file

      USE iotk_module
      USE io_files, ONLY: iuni_nml, nmlfile
      IMPLICIT NONE
      CHARACTER(LEN=256) :: attr

      IF( mpime == root ) THEN

          OPEN( UNIT=iuni_nml, FILE=nmlfile, STATUS='UNKNOWN')
          WRITE(iuni_nml, * ) '<?xml version="1.0" encoding="UTF-8"?>'

          CALL iotk_write_begin( iuni_nml, "input" )

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
            CALL iotk_write_dat( iuni_nml, "imr", imr )
            CALL iotk_write_dat( iuni_nml, "tpr", tpr )
            CALL iotk_write_dat( iuni_nml, "tdump", tdump )
            CALL iotk_write_dat( iuni_nml, "nfil", nfil )
            CALL iotk_write_dat( iuni_nml, "max_seconds", max_seconds )
            CALL iotk_write_dat( iuni_nml, "tau", tau )
            CALL iotk_write_dat( iuni_nml, "formatted_output", formatted_output )
          CALL iotk_write_end( iuni_nml, "control" )

          CALL iotk_write_begin( iuni_nml, "model" )
            CALL iotk_write_dat( iuni_nml, "icpc", icpc )
            CALL iotk_write_dat( iuni_nml, "irex", irex )
            CALL iotk_write_dat( iuni_nml, "gas_viscosity", gas_viscosity )
            CALL iotk_write_dat( iuni_nml, "part_viscosity", part_viscosity )
            CALL iotk_write_dat( iuni_nml, "iss", iss )
            CALL iotk_write_dat( iuni_nml, "repulsive_model",repulsive_model )
            CALL iotk_write_dat( iuni_nml, "iturb", iturb )
            CALL iotk_write_dat( iuni_nml, "modturbo", modturbo )
            CALL iotk_write_dat( iuni_nml, "cmut", cmut )
            CALL iotk_write_dat( iuni_nml, "gravx", gravx )
            CALL iotk_write_dat( iuni_nml, "gravy", gravy )
            CALL iotk_write_dat( iuni_nml, "gravz", gravz )
            CALL iotk_write_dat( iuni_nml, "ngas", ngas )
            CALL iotk_write_dat( iuni_nml, "density_specified", density_specified )
          CALL iotk_write_end( iuni_nml, "model" )


          CALL iotk_write_begin( iuni_nml, "mesh" )
            CALL iotk_write_dat( iuni_nml, "nx", nx )
            CALL iotk_write_dat( iuni_nml, "ny", ny )
            CALL iotk_write_dat( iuni_nml, "nz", nz )
            CALL iotk_write_dat( iuni_nml, "itc", itc )
            CALL iotk_write_dat( iuni_nml, "iuni", iuni )
            CALL iotk_write_dat( iuni_nml, "dx0", dx0 )
            CALL iotk_write_dat( iuni_nml, "dy0", dy0 )
            CALL iotk_write_dat( iuni_nml, "dz0", dz0 )
            CALL iotk_write_dat( iuni_nml, "zzero", zzero )
            CALL iotk_write_dat( iuni_nml, "center_x", center_x )
            CALL iotk_write_dat( iuni_nml, "center_y", center_y )
            CALL iotk_write_dat( iuni_nml, "alpha_x", alpha_x )
            CALL iotk_write_dat( iuni_nml, "alpha_y", alpha_y )
            CALL iotk_write_dat( iuni_nml, "alpha_z", alpha_z )
            CALL iotk_write_dat( iuni_nml, "dxmin", dxmin )
            CALL iotk_write_dat( iuni_nml, "dymin", dymin )
            CALL iotk_write_dat( iuni_nml, "dzmin", dzmin )
            CALL iotk_write_dat( iuni_nml, "dxmax", dxmax )
            CALL iotk_write_dat( iuni_nml, "dymax", dymax )
            CALL iotk_write_dat( iuni_nml, "dzmax", dzmax )
            CALL iotk_write_dat( iuni_nml, "n0x", n0x )
            CALL iotk_write_dat( iuni_nml, "n0y", n0y )
            CALL iotk_write_dat( iuni_nml, "n0z", n0z )
            CALL iotk_write_dat( iuni_nml, "domain_x", domain_x )
            CALL iotk_write_dat( iuni_nml, "domain_y", domain_y )
            CALL iotk_write_dat( iuni_nml, "domain_z", domain_z )
            CALL iotk_write_dat( iuni_nml, "maxbeta", maxbeta )
            CALL iotk_write_dat( iuni_nml, "grigen", grigen )
            CALL iotk_write_dat( iuni_nml, "mesh_partition", mesh_partition )
          CALL iotk_write_end( iuni_nml, "mesh" )

          CALL iotk_write_begin( iuni_nml, "boundaries" )
            CALL iotk_write_dat( iuni_nml, "east", east )
            CALL iotk_write_dat( iuni_nml, "west", west )
            CALL iotk_write_dat( iuni_nml, "south", south )
            CALL iotk_write_dat( iuni_nml, "north", north )
            CALL iotk_write_dat( iuni_nml, "top", top )
            CALL iotk_write_dat( iuni_nml, "bottom", bottom )
            CALL iotk_write_dat( iuni_nml, "immb", immb )
            CALL iotk_write_dat( iuni_nml, "ibl", ibl )
          CALL iotk_write_end( iuni_nml, "boundaries" )

          CALL iotk_write_begin( iuni_nml, "topography" )
            CALL iotk_write_begin( iuni_nml, "dem_file" )
              WRITE( iuni_nml, * ) dem_file
            CALL iotk_write_end( iuni_nml, "dem_file" )
            CALL iotk_write_dat( iuni_nml, "itp", itp )
            CALL iotk_write_dat( iuni_nml, "ismt", ismt )
            CALL iotk_write_dat( iuni_nml, "iavv", iavv )
            CALL iotk_write_dat( iuni_nml, "nocrater", nocrater )
            CALL iotk_write_dat( iuni_nml, "itrans", itrans )
            CALL iotk_write_dat( iuni_nml, "rim_quota", rim_quota )
            CALL iotk_write_dat( iuni_nml, "filtersize", filtersize )
            CALL iotk_write_dat( iuni_nml, "cellsize", cellsize )
            CALL iotk_write_dat( iuni_nml, "zrough", zrough )
          CALL iotk_write_end( iuni_nml, "topography" )

          CALL iotk_write_begin( iuni_nml, "inlet" )
            CALL iotk_write_dat( iuni_nml, "ivent", ivent )
            CALL iotk_write_dat( iuni_nml, "iali", iali )
            CALL iotk_write_dat( iuni_nml, "ipro", ipro )
            CALL iotk_write_begin( iuni_nml, "rad_file" )
              WRITE( iuni_nml, * ) rad_file
            CALL iotk_write_end( iuni_nml, "rad_file" )
            CALL iotk_write_dat( iuni_nml, "wrat", wrat )
            CALL iotk_write_dat( iuni_nml, "irand", irand )
            CALL iotk_write_dat( iuni_nml, "xvent", xvent )
            CALL iotk_write_dat( iuni_nml, "yvent", yvent )
            CALL iotk_write_dat( iuni_nml, "vent_radius", vent_radius )
            CALL iotk_write_dat( iuni_nml, "base_radius", base_radius )
            CALL iotk_write_dat( iuni_nml, "crater_radius", crater_radius )
            CALL iotk_write_dat( iuni_nml, "u_gas", u_gas )
            CALL iotk_write_dat( iuni_nml, "v_gas", v_gas )
            CALL iotk_write_dat( iuni_nml, "w_gas", w_gas )
            CALL iotk_write_dat( iuni_nml, "p_gas", p_gas )
            CALL iotk_write_dat( iuni_nml, "t_gas", t_gas )
            CALL iotk_write_dat( iuni_nml, "u_solid", u_solid )
            CALL iotk_write_dat( iuni_nml, "v_solid", v_solid )
            CALL iotk_write_dat( iuni_nml, "w_solid", w_solid )
            CALL iotk_write_dat( iuni_nml, "ep_solid", ep_solid )
            CALL iotk_write_dat( iuni_nml, "t_solid", t_solid )
            CALL iotk_write_dat( iuni_nml, "vent_O2", vent_O2 )
            CALL iotk_write_dat( iuni_nml, "vent_N2", vent_N2 )
            CALL iotk_write_dat( iuni_nml, "vent_CO2", vent_CO2 )
            CALL iotk_write_dat( iuni_nml, "vent_H2", vent_H2 )
            CALL iotk_write_dat( iuni_nml, "vent_H2O", vent_H2O )
            CALL iotk_write_dat( iuni_nml, "vent_Air", vent_Air )
            CALL iotk_write_dat( iuni_nml, "vent_SO2", vent_SO2 )
          CALL iotk_write_end( iuni_nml, "inlet" )

          CALL iotk_write_begin( iuni_nml, "dome" )
            CALL iotk_write_dat( iuni_nml, "idome", idome )
            CALL iotk_write_dat( iuni_nml, "idw", idw )
            CALL iotk_write_dat( iuni_nml, "xdome", xdome )
            CALL iotk_write_dat( iuni_nml, "ydome", ydome )
            CALL iotk_write_dat( iuni_nml, "dome_volume", dome_volume )
            CALL iotk_write_dat( iuni_nml, "temperature", temperature )
            CALL iotk_write_dat( iuni_nml, "overpressure", overpressure )
            CALL iotk_write_dat( iuni_nml, "particle_fraction", particle_fraction )
            CALL iotk_write_dat( iuni_nml, "gas_flux", gas_flux )
            CALL iotk_write_dat( iuni_nml, "permeability", permeability )
            CALL iotk_write_dat( iuni_nml, "dome_gasvisc", dome_gasvisc )
            CALL iotk_write_dat( iuni_nml, "dome_O2", dome_O2 )
            CALL iotk_write_dat( iuni_nml, "dome_N2", dome_N2 )
            CALL iotk_write_dat( iuni_nml, "dome_CO2", dome_CO2 )
            CALL iotk_write_dat( iuni_nml, "dome_H2", dome_H2 )
            CALL iotk_write_dat( iuni_nml, "dome_H2O", dome_H2O )
            CALL iotk_write_dat( iuni_nml, "dome_Air", dome_Air )
            CALL iotk_write_dat( iuni_nml, "dome_SO2", dome_SO2 )
          CALL iotk_write_end( iuni_nml, "dome" )

          CALL iotk_write_begin( iuni_nml, "atmosphere" )
            CALL iotk_write_dat( iuni_nml, "wind_x", wind_x )
            CALL iotk_write_dat( iuni_nml, "wind_y", wind_y )
            CALL iotk_write_dat( iuni_nml, "wind_z", wind_z )
            CALL iotk_write_dat( iuni_nml, "stratification", stratification )
            CALL iotk_write_dat( iuni_nml, "p_ground", p_ground )
            CALL iotk_write_dat( iuni_nml, "t_ground", t_ground )
            CALL iotk_write_dat( iuni_nml, "void_fraction", void_fraction )
            CALL iotk_write_dat( iuni_nml, "max_packing", max_packing )
            CALL iotk_write_dat( iuni_nml, "atm_O2", atm_O2 )
            CALL iotk_write_dat( iuni_nml, "atm_N2", atm_N2 )
            CALL iotk_write_dat( iuni_nml, "atm_CO2", atm_CO2 )
            CALL iotk_write_dat( iuni_nml, "atm_H2", atm_H2 )
            CALL iotk_write_dat( iuni_nml, "atm_H2O", atm_H2O )
            CALL iotk_write_dat( iuni_nml, "atm_Air", atm_Air )
            CALL iotk_write_dat( iuni_nml, "atm_SO2", atm_SO2 )
            CALL iotk_write_dat( iuni_nml, "troposphere_z",troposphere_z)
            CALL iotk_write_dat( iuni_nml, "tropopause_z",tropopause_z)
            CALL iotk_write_dat( iuni_nml, "lower_stratosphere_z",lower_stratosphere_z)
            CALL iotk_write_dat( iuni_nml, "upper_stratosphere_z",upper_stratosphere_z)
            CALL iotk_write_dat( iuni_nml, "ozone_layer_z",ozone_layer_z)
            CALL iotk_write_dat( iuni_nml, "lower_mesosphere_z",lower_mesosphere_z)
            CALL iotk_write_dat( iuni_nml, "upper_mesosphere_z",upper_mesosphere_z)
            CALL iotk_write_dat( iuni_nml, "troposphere_grad",troposphere_grad)
            CALL iotk_write_dat( iuni_nml, "tropopause_grad",tropopause_grad)
            CALL iotk_write_dat( iuni_nml, "lower_stratosphere_grad",lower_stratosphere_grad)
            CALL iotk_write_dat( iuni_nml, "upper_stratosphere_grad",upper_stratosphere_grad)
            CALL iotk_write_dat( iuni_nml, "ozone_layer_grad",ozone_layer_grad)
            CALL iotk_write_dat( iuni_nml, "lower_mesosphere_grad",lower_mesosphere_grad)
            CALL iotk_write_dat( iuni_nml, "upper_mesosphere_grad",upper_mesosphere_grad)
          CALL iotk_write_end( iuni_nml, "atmosphere" )

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

          CALL iotk_write_begin( iuni_nml, "numeric" )
            CALL iotk_write_dat( iuni_nml, "rungekut", rungekut )
            CALL iotk_write_dat( iuni_nml, "beta", beta )
            CALL iotk_write_dat( iuni_nml, "muscl", muscl )
            CALL iotk_write_dat( iuni_nml, "mass_limiter", mass_limiter )
            CALL iotk_write_dat( iuni_nml, "vel_limiter", vel_limiter )
            CALL iotk_write_dat( iuni_nml, "inmax", inmax )
            CALL iotk_write_dat( iuni_nml, "delg", delg )
            CALL iotk_write_dat( iuni_nml, "lim_type", lim_type )
            CALL iotk_write_dat( iuni_nml, "maxout", maxout )
            CALL iotk_write_dat( iuni_nml, "optimization", optimization )
            CALL iotk_write_dat( iuni_nml, "omega", omega )
            CALL iotk_write_dat( iuni_nml, "tforce", tforce)
            CALL iotk_write_dat( iuni_nml, "vforce", vforce)
            CALL iotk_write_dat( iuni_nml, "implicit_fluxes",   &
                                          & implicit_fluxes )
            CALL iotk_write_dat( iuni_nml, "implicit_enthalpy", &
                                          & implicit_enthalpy )
            CALL iotk_write_dat( iuni_nml, "update_eosg", &
                                          & update_eosg )
            CALL iotk_write_dat( iuni_nml, "rlim", rlim )
            CALL iotk_write_dat( iuni_nml, "flim", flim )
          CALL iotk_write_end( iuni_nml, "numeric" )

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

          CALL iotk_write_end( iuni_nml, "input" )

          CLOSE( UNIT=iuni_nml )

      END IF

      RETURN
      END SUBROUTINE write_xml
!----------------------------------------------------------------------
      END SUBROUTINE input
!----------------------------------------------------------------------
      SUBROUTINE initc

      USE atmospheric_conditions, ONLY: atm_ygc, layer
      USE control_flags, ONLY: job_type
      USE dimensions
      USE dome_conditions, ONLY: dome_ygc, idome
      USE flux_limiters, ONLY: lv, lm
      USE grid, ONLY: dx, dy, dz, itc
      USE grid, ONLY: iob, zzero, grigen
      USE immersed_boundaries, ONLY: immb
      USE initial_conditions, ONLY: epsob, tpob, ygcob,   &
     &     ugob, vgob, wgob, upob, vpob, wpob, pob, tgob, epob
      USE particles_constants, ONLY: rl, inrl, kap, cmus, phis, cps, dk
      USE vent_conditions, ONLY: vent_ygc, ivent
      USE volcano_topography, ONLY: itp

      IMPLICIT NONE

      INTEGER :: ig, is
!
! ... flags compatibility (immersed boundaries are used only with 
! ... implicit topography
!
      IF (itp < 1) immb = 0
      IF (idome >= 1) ivent = 0
!
! ... numeric
!
      lv = vel_limiter
      lm = mass_limiter
!
! ... mesh
!
      IF( grigen == 0) THEN
        dx(1:nx) = delta_x(1:nx)
        IF( job_type == '3D' ) THEN
          dy(1:ny) = delta_y(1:ny)
        ELSE IF( job_type == '2D' ) THEN
          dy(1:ny) = 0.D0
        END IF
        dz(1:nz) = delta_z(1:nz)
      END IF
!
! ... specified flows
!
      IF (no > 0) THEN
        !
        ! ... blocks specification
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
        !
        ! ... specified flow
        !
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
      END IF
!
! ... Atmospheric composition 
!
      atm_ygc(1) = atm_O2
      atm_ygc(2) = atm_N2
      atm_ygc(3) = atm_CO2
      atm_ygc(4) = atm_H2
      atm_ygc(5) = atm_H2O
      atm_ygc(6) = atm_Air
      atm_ygc(7) = atm_SO2
!
! ... Initialize atmospheric layers
!
      layer(1)%name = 'Troposphere'   
      layer(2)%name = 'Tropopause'    
      layer(3)%name = 'Lower_Stratosphere' 
      layer(4)%name = 'Upper_Stratosphere' 
      layer(5)%name = 'Ozone_layer'   
      layer(6)%name = 'Lower_Mesosphere'   
      layer(7)%name = 'Upper_Mesosphere'   
!
! ... Top of the layer (height a.s.l.)
!
      layer(1)%ztop = troposphere_z
      layer(2)%ztop = tropopause_z
      layer(3)%ztop = lower_stratosphere_z
      layer(4)%ztop = upper_stratosphere_z
      layer(5)%ztop = ozone_layer_z
      layer(6)%ztop = lower_mesosphere_z
      layer(7)%ztop = upper_mesosphere_z
!
! ... Temperature gradient (T is assumed to vary linearly)
!
      layer(1)%gradt = troposphere_grad
      layer(2)%gradt = tropopause_grad
      layer(3)%gradt = lower_stratosphere_grad
      layer(4)%gradt = upper_stratosphere_grad
      layer(5)%gradt = ozone_layer_grad
      layer(6)%gradt = lower_mesosphere_grad
      layer(7)%gradt = upper_mesosphere_grad
!
! ... Gas components at vent
!
      IF (ivent >= 1) THEN
        vent_ygc(1) = vent_O2
        vent_ygc(2) = vent_N2
        vent_ygc(3) = vent_CO2
        vent_ygc(4) = vent_H2
        vent_ygc(5) = vent_H2O
        vent_ygc(6) = vent_Air
        vent_ygc(7) = vent_SO2
      ELSE
        vent_ygc = 0.D0
      END IF
!
! ... Gas components at dome
!
      IF (idome >= 1) THEN
              dome_ygc(1) = dome_O2
              dome_ygc(2) = dome_N2
              dome_ygc(3) = dome_CO2
              dome_ygc(4) = dome_H2
              dome_ygc(5) = dome_H2O
              dome_ygc(6) = dome_Air
              dome_ygc(7) = dome_SO2
      ELSE
              dome_ygc = 0.D0
      END IF
!
! ... particle properties
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
