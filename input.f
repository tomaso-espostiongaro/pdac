!-----------------------------------------------------------------------------------------
   MODULE input_module
!-----------------------------------------------------------------------------------------

      USE dimensions, ONLY: max_nsolid, ngas, nroughx, max_size, max_nblock, max_ngas, nr, nz

      REAL*8 :: diameter(max_nsolid)
      REAL*8 :: density(max_nsolid)
      REAL*8 :: sphericity(max_nsolid)
      REAL*8 :: viscosity(max_nsolid)
      REAL*8 :: specific_heat(max_nsolid)
      REAL*8 :: thermal_conductivity(max_nsolid)
      REAL*8 :: dr0, dz0
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
      INTEGER :: block_bounds(4,max_nblock)
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

      USE atmosphere, ONLY: v0, u0, p0, temp0, uk0, vk0, ep0, epsmx0, gravx, gravz
      USE eulerian_flux, ONLY: beta, muscl
      USE gas_constants, ONLY: phij, ckg, mmug, mmugs, mmugek, gmw
      USE grid, ONLY: dz, dr, itc
      USE grid, ONLY: nso, iob
      USE iterative_solver, ONLY: inmax, maxout, omega
      USE output_dump, ONLY: nfil, outp
      USE parallel, ONLY: mpime, root
      USE particles_constants, ONLY: rl, inrl, kap, &
     &     cmus, phis, cps, dk, nsolid
      USE phases_matrix, ONLY: rlim
      USE reactions, ONLY: irex
      USE roughness_module, ONLY: zrough, allocate_roughness, roughness
      USE initial_conditions, ONLY: setup, epsob, vpob, tpob, ygc0, &
     &    ygcob, upob, vgob, ugob, pob, tgob, epob, lpr, zzero, &
     &    bounds_setup
      USE time_parameters, ONLY: time, tstop, dt, tpr, tdump, itd, & 
     &                            timestart, rungekut
      USE turbulence, ONLY: iturb, cmut, iss, modturbo
!
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: iunit

      NAMELIST / control / run_name, restart_mode, time, tstop, dt, lpr, tpr, &
        tdump, nfil, irex, iss, iturb, modturbo, cmut, rlim, gravx, gravz, &
        ngas
!
      NAMELIST / mesh / nr, nz, itc, iuni, dr0, dz0
!
      NAMELIST / particles / nsolid, diameter, density, sphericity, &
        viscosity, specific_heat, thermal_conductivity
!
      NAMELIST / numeric / rungekut, beta, muscl, inmax, maxout, omega
!
      INTEGER :: i, j, k, n, m, kg
      CHARACTER(LEN=80) :: card
      LOGICAL :: tend
!
!    Sets default values

! ... Control

      run_name = 'run2d'
      restart_mode = 'from_scratch'
      time = 0.0d0
      tstop = 1.0d0
      dt = 0.01d0
      lpr = 1
      tpr = 1.0d0
      tdump = 20.0d0
      nfil = 0
      irex = 1
      iss  = 0
      iturb = 1
      modturbo = 1
      cmut = 0.1d0
      rlim = 1.0d-8
      gravx = 0.0d0
      gravz = -981.0d0
      ngas = 7

! ... Mesh

      nr = 100
      nz = 100
      itc = 0
      iuni = 0
      dr0  = 1000.0d0
      dz0  = 1000.0d0
      zzero = 0.0d0

! ... Particles
 
      nsolid = 2
      diameter = 0.01
      density = 2.7
      sphericity = 1.0
      viscosity = 5.0
      specific_heat = 1.2d7
      thermal_conductivity = 2.0d5

! ... Numeric

      rungekut = 1
      beta = 0.25
      muscl = 0.0
      inmax = 8
      maxout = 5000
      omega = 1.0

!
! reading of input file
!

      IF(mpime .EQ. root) THEN
        READ(iunit, control) 
      END IF
!
      CALL bcast_character(run_name,80,root)
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

      SELECT CASE ( TRIM(restart_mode) )
        CASE ('from_scratch')
          itd = 1 
        CASE ('restart', 'default' )
          itd = 2
        CASE DEFAULT
          CALL error(' input ',' unknown restart_mode '//TRIM(restart_mode), 1 )
      END SELECT

      IF(mpime .EQ. root) THEN
        READ(iunit, mesh) 
      END IF

      CALL bcast_integer(nr,1,root)
      CALL bcast_integer(nz,1,root)
      CALL bcast_integer(itc,1,root)
      CALL bcast_integer(iuni,1,root)
      CALL bcast_real(dr0,1,root)
      CALL bcast_real(dz0,1,root)
      CALL bcast_real(zzero,1,root)

      IF(mpime .EQ. root) THEN
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

      IF(mpime .EQ. root) THEN
        READ(iunit, numeric) 
      END IF

      CALL bcast_integer(rungekut,1,root)
      CALL bcast_real(beta,1,root)
      CALL bcast_real(muscl,1,root)
      CALL bcast_integer(inmax,1,root)
      CALL bcast_integer(maxout,1,root)
      CALL bcast_real(omega,1,root)

      CALL allocate_roughness( zrough, nroughx )
      tend = .FALSE.
      IF(mpime .EQ. root) THEN
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
      IF(mpime .EQ. root) THEN
        mesh_search: DO
          READ(5,*,END=200) card
          IF( TRIM(card) == 'MESH' ) THEN
            EXIT mesh_search
          END IF
        END DO mesh_search

        READ(5,*) (delta_r(i),i=1,nr)
        READ(5,*) (delta_z(j),j=1,nz)
        IF (iuni == 0) THEN
          delta_r = dr0
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
      CALL bcast_real(delta_z,nz,root)


      tend = .FALSE.
      IF(mpime .EQ. root) THEN
        fixed_flows_search: DO
          READ(5,*,END=300) card
          IF( TRIM(card) == 'FIXED_FLOWS' ) THEN
            EXIT fixed_flows_search
          END IF
        END DO fixed_flows_search

        READ(5,*) number_of_block
        DO n = 1, number_of_block
          READ(5,*) block_type(n),(block_bounds(m,n),m=1,4)
          IF( block_type(n) == 1 .OR. block_type(n) == 5) THEN
            READ(5,*) fixed_vgas_r(n), fixed_vgas_z(n), fixed_pressure(n), fixed_gaseps(n), fixed_gastemp(n)
            READ(5,*) ( fixed_vpart_r(k,n), fixed_vpart_z(k,n), fixed_parteps(k,n), fixed_parttemp(k,n), k=1, nsolid)
            READ(5,*) ( fixed_gasconc(kg,n), kg=1, ngas )
          ENDIF
        END DO

        GOTO 310
 300    tend = .TRUE.
 310    continue
      END IF
!
      CALL bcast_logical(tend, 1, root)
      IF( tend ) THEN
        CALL error( ' input ', ' MESH card not found ', 1 )
      END IF
      CALL bcast_integer(number_of_block, 1, root)
      CALL bcast_integer(block_type, SIZE(block_type), root)
      CALL bcast_integer(block_bounds, SIZE(block_bounds), root)
      CALL bcast_real(fixed_vgas_r, SIZE(fixed_vgas_r), root)
      CALL bcast_real(fixed_vgas_z, SIZE(fixed_vgas_z), root)
      CALL bcast_real(fixed_pressure, SIZE(fixed_pressure), root)
      CALL bcast_real(fixed_gaseps, SIZE(fixed_gaseps), root)
      CALL bcast_real(fixed_gastemp, SIZE(fixed_gastemp), root)
      CALL bcast_real(fixed_vpart_r, SIZE(fixed_vpart_r), root)
      CALL bcast_real(fixed_vpart_z, SIZE(fixed_vpart_z), root)
      CALL bcast_real(fixed_parteps, SIZE(fixed_parteps), root)
      CALL bcast_real(fixed_parttemp, SIZE(fixed_parttemp), root)
      CALL bcast_real(fixed_gasconc, SIZE(fixed_gasconc), root)


      tend = .FALSE.
      IF(mpime .EQ. root) THEN
        initial_conditions_search: DO
          READ(5,*,END=400) card
          IF( TRIM(card) == 'INITIAL_CONDITIONS' ) THEN
            EXIT initial_conditions_search
          END IF
        END DO initial_conditions_search

        READ(5,*) initial_vgas_r, initial_vgas_z, initial_pressure, initial_void_fraction, max_packing, initial_temperature
        READ(5,*) initial_vpart_r, initial_vpart_z
        READ(5,*) (initial_gasconc(kg), kg=1, ngas)

        GOTO 410
 400    tend = .TRUE.
 410    continue
      END IF
!
      CALL bcast_logical(tend, 1, root)
      IF( tend ) THEN
        CALL error( ' input ', ' MESH card not found ', 1 )
      END IF
      CALL bcast_real(initial_vgas_r,1,root)
      CALL bcast_real(initial_vgas_z,1,root)
      CALL bcast_real(initial_pressure,1,root)
      CALL bcast_real(initial_void_fraction,1,root)
      CALL bcast_real(max_packing,1,root)
      CALL bcast_real(initial_temperature,1,root)
      CALL bcast_real(initial_vpart_r,1,root)
      CALL bcast_real(initial_vpart_z,1,root)
      CALL bcast_real(initial_gasconc, SIZE(initial_gasconc),root)

      END SUBROUTINE

  END MODULE input_module
