   MODULE input_module

      USE dimensions, ONLY: max_nsolid, ngas

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

     CONTAINS
!
      SUBROUTINE input( iunit )

      USE atmosphere, ONLY: v0, u0, p0, temp0, uk0, vk0, ep0, epsmx0, gravx, gravz
      USE gas_constants, ONLY: phij, ckg, mmug, mmugs, mmugek, gmw
      USE gas_solid_viscosity, ONLY: icoh
      USE grid, ONLY: dz, dr, itc, ib2, ib1, ib, jb2, jb1, jb
      USE grid, ONLY: no, nso, iob
      USE iterative_solver, ONLY: inmax, maxout, omega
      USE output_dump, ONLY: nfil, outp
      USE parallel, ONLY: mpime, root
      USE particles_constants, ONLY: rl, inrl, kap, &
     &     cmus, phis, cps, dk, nsolid
      USE phases_matrix, ONLY: rlim
      USE reactions, ONLY: irex
      USE roughness, ONLY: roucha, zrough, ir
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
        tdump, nfil, icoh, irex, iss, iturb, modturbo, cmut, rlim, gravx, gravz, &
        ngas
!
      NAMELIST / mesh / ib2, jb2, itc, iuni, dr0, dz0
!
      NAMELIST / particles / nsolid, diameter, density, sphericity, &
        viscosity, specific_heat, thermal_conductivity
!
      NAMELIST / numeric / rungekut, inmax, maxout, omega
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
      icoh = 0
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

      ib2 = 100
      jb2 = 100
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
      CALL bcast_integer(icoh,1,root)
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

      CALL bcast_integer(ib2,1,root)
      CALL bcast_integer(jb2,1,root)
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
      CALL bcast_integer(inmax,1,root)
      CALL bcast_integer(maxout,1,root)
      CALL bcast_real(omega,1,root)

      END SUBROUTINE

  END MODULE input_module
