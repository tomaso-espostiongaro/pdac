!******************************************************************************
! PROGRAM: pdac (pyroclastic dispersion analysis code)
! description: parallel multiphase and multicomponent fluid dynamic code
!              for the simulation of pyroclastic dispersion
!              processes
! version: 3.0
! date: June 2003
! authors: A.Neri, G.Macedonio, D.Gidaspow, T. Esposti Ongaro
! parallelized by: T. Esposti Ongaro, C. Cavazzoni, A. Neri
! 3d version implemented by: T. Esposti Ongaro, C. Cavazzoni
!******************************************************************************
!
      PROGRAM pdac

      USE atmosphere, ONLY: v0, u0, w0, p0, temp0, us0, vs0, ws0, &
     &                      ep0, epsmx0, gravx, gravz
      USE dimensions
      USE eos_gas, ONLY: allocate_eosg
      USE gas_constants, ONLY: allocate_gas_constants, present_gas
      USE gas_solid_density, ONLY: allocate_density
      USE gas_solid_velocity, ONLY: allocate_velocity
      USE gas_solid_temperature, ONLY: allocate_temperature
      USE gas_solid_viscosity, ONLY: allocate_viscosity
      USE grid, ONLY: dx, dy, dz, dr, itc
      USE grid, ONLY: iob, flic, partition, ghost, &
     &    allocate_blbody, allocate_grid
      USE io_restart, ONLY: taperd, tapewr
      USE iterative_solver, ONLY: inmax, maxout, omega
      USE output_dump, ONLY: nfil, recover_2d
      USE parallel, ONLY: parallel_startup, parallel_hangup, &
     &    mpime, root
      USE particles_constants, ONLY: rl, inrl, kap, &
     &     cmus, phis, cps, dk, allocate_part_constants
      USE phases_matrix, ONLY: rlim, allocate_matrix
      USE pressure_epsilon, ONLY: allocate_press_eps
      USE roughness_module, ONLY: zrough, deallocate_roughness
      USE initial_conditions, ONLY: setup, epsob, tpob, ygc0, ygcob, &
     &     ugob, vgob, wgob, upob, vpob, wpob, pob, tgob, epob, lpr, & 
     &     zzero, allocate_setup
      USE specific_heat_module, ONLY: allocate_hcapgs
      USE time_parameters, ONLY: time, tstop, dt, tpr, tdump, itd, & 
     &                            timestart, rungekut
      USE turbulence_model, ONLY: allocate_turbo
      USE environment, ONLY: cpclock, timing
      USE input_module
      USE control_flags, ONLY: job_type
!
      IMPLICIT NONE
      CHARACTER(LEN=11) :: errnb, testnb
      CHARACTER(LEN=8) :: inputfile, logfile, errorfile, testfile
      CHARACTER(LEN=3) :: procnum
!
      INTEGER :: inputunit, logunit, errorunit, testunit
      INTEGER :: ig
      REAL*8 :: s0, s1, s2, s3, s4
      REAL*8 :: timtot, timprog, timdist, timsetup, timinit
      LOGICAL :: debug = .TRUE.
!
      IF(timing) s0 = cpclock()

!
! ... initialize parallel environment

      CALL parallel_startup  

! ... Initialize the IBM HW performance monitor
!      call f_hpminit( mpime, 'pdac' )
!
! ... I/O files
!
      inputunit = 5
      logunit   = 6
      testunit  = 7
      errorunit = 8
      inputfile = 'pdac.dat'
      logfile = 'pdac.log'
      testfile = 'pdac.tst'
      errorfile = 'pdac.err'
      errnb = errorfile//procnum(mpime)
      testnb = testfile//procnum(mpime)
      IF(mpime .EQ. root) THEN
        OPEN(UNIT=inputunit, FILE=inputfile, STATUS='UNKNOWN')
        IF( .NOT. debug ) OPEN(UNIT=logunit,   FILE=logfile,   STATUS='UNKNOWN')
        OPEN(UNIT=testunit,  FILE=testfile,  STATUS='UNKNOWN')
        OPEN(UNIT=errorunit, FILE=errorfile, STATUS='UNKNOWN')
      ELSE
        OPEN(UNIT=testunit,  FILE=testnb,    STATUS='UNKNOWN')
        OPEN(UNIT=errorunit, FILE=errnb,     STATUS='UNKNOWN')
      END IF
!
! ... Read Input file
!
      CALL input(inputunit)

! ... allocation of arrays ...(dependence on nr, nx,ny,nz)
      CALL allocate_grid
!
      IF( job_type == '2D' ) THEN
        dr(1:nr) = delta_r(1:nr)
      ELSE IF( job_type == '3D' ) THEN
        dx(1:nx) = delta_x(1:nx)
        dy(1:ny) = delta_y(1:ny)
      END IF
      dz(1:nz) = delta_z(1:nz)

      no = number_of_block
!
! ... set dimensions ...
      nphase=nsolid+1

! ... allocation of arrays ...(dependence on no, ngas, nsolid, nr, nz)
      CALL allocate_blbody
      CALL allocate_gas_constants
      CALL allocate_matrix
      CALL allocate_part_constants
      CALL allocate_setup

      iob(1:no)%typ = block_type(1:no)

      IF( job_type == '2D' ) THEN
        iob(1:no)%rlo = block_bounds(1,1:no)
        iob(1:no)%rhi = block_bounds(2,1:no)
      ELSE IF ( job_type == '3D' ) THEN
        iob(1:no)%xlo = block_bounds(1,1:no)
        iob(1:no)%xhi = block_bounds(2,1:no)
        iob(1:no)%ylo = block_bounds(3,1:no)
        iob(1:no)%yhi = block_bounds(4,1:no)
      END IF
      iob(1:no)%zlo = block_bounds(5,1:no)
      iob(1:no)%zhi = block_bounds(6,1:no)

      IF( job_type == '2D' ) THEN
        ugob(1:no)  = fixed_vgas_r(1:no)
      ELSE
        ugob(1:no)  = fixed_vgas_x(1:no)
        vgob(1:no)  = fixed_vgas_y(1:no)
      END IF
      wgob(1:no)  = fixed_vgas_z(1:no)
      pob(1:no)  = fixed_pressure(1:no)
      epob(1:no)  = fixed_gaseps(1:no)
      tgob(1:no)  = fixed_gastemp(1:no)
      IF( job_type == '2D' ) THEN
        upob(1:nsolid,1:no) = fixed_vpart_r(1:nsolid,1:no)
      ELSE
        upob(1:nsolid,1:no) = fixed_vpart_x(1:nsolid,1:no)
        vpob(1:nsolid,1:no) = fixed_vpart_y(1:nsolid,1:no)
      END IF
      wpob(1:nsolid,1:no) = fixed_vpart_z(1:nsolid,1:no)
      epsob(1:nsolid,1:no) = fixed_parteps(1:nsolid,1:no)
      tpob(1:nsolid,1:no) = fixed_parttemp(1:nsolid,1:no)
      ygcob(1:ngas,1:no) = fixed_gasconc(1:ngas,1:no)

      IF( job_type == '2D' ) THEN
        u0 = initial_vgas_r
      ELSE
        u0 = initial_vgas_x
        v0 = initial_vgas_y
      END IF
      w0 = initial_vgas_z

      p0 = initial_pressure
      ep0 = initial_void_fraction
      epsmx0 = max_packing
      temp0 = initial_temperature

      IF( job_type == '2D' ) THEN
        us0 = initial_vpart_r
      ELSE
        us0 = initial_vpart_x
        vs0 = initial_vpart_y
      END IF
      ws0 = initial_vpart_z

      ygc0(1:ngas) = initial_gasconc(1:ngas)
!
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

! ... set start time
      timestart = time
!
! ... Set cell-type flags
!
      CALL flic
      IF(timing) s1 = cpclock()
!
! ... Domain decomposition for parallelization 
!
      CALL partition
      IF(timing) s2 = cpclock()
!
! ... Setting ghost cells for parallel data exchange
! ... and the indexes
!
      CALL ghost
      IF(timing) s3 = cpclock()
!
      CALL allocate_velocity
      CALL allocate_density
      CALL allocate_press_eps
      CALL allocate_temperature
      CALL allocate_eosg
      CALL allocate_viscosity
      CALL allocate_hcapgs
      CALL allocate_turbo

!
! ... Read restart file
!

      IF(itd == 2) THEN 
        CALL taperd
      ELSE IF (itd >= 3) THEN
        CALL recover_2d
      END IF

!
! ... Set initial conditions
!
      CALL setup
!
! ... Time advancement loop
!
      CALL prog
!
      IF(timing ) THEN
        s4 = cpclock()
        timtot   = (s4 - s0)/1000.D0
        timprog  = (s4 - s3)/1000.D0
        timdist  = (s3 - s2)/1000.D0
        timsetup = (s2 - s1)/1000.D0
        timinit  = (s1 - s0)/1000.D0
        WRITE(7,900)
        WRITE(7,999) timinit,timsetup,timdist,timprog,timtot
999     FORMAT(6(1X,F9.3))
900     FORMAT('  Init      Part ',  '    Ghost      Prog      Total')
      END IF

! ... terminate the IBM HW performance monitor session
!      call f_hpmterminate( mpime )
!
      CLOSE(5)
      IF( .NOT. debug ) CLOSE(6)
      CLOSE(7)
      CLOSE(8)
!
      CALL deallocate_roughness( zrough )
! 
! ... Finalize parallel environment
!
      CALL parallel_hangup
!
      STOP
!
!**************************************************************
      END PROGRAM pdac
!**************************************************************
