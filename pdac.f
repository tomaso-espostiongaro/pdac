!******************************************************************************
! PROGRAM: pdac (pyroclastic dispersion analysis code)
! description: parallel multiphase and multicomponent fluid dynamic code
!              for the simulation of pyroclastic dispersion
!              processes
! version: 3.0
! date: Jan 2005
! authors: A.Neri, G.Macedonio, D.Gidaspow, T. Esposti Ongaro
! parallelized by: T. Esposti Ongaro, C. Cavazzoni, A. Neri
! 3d version implemented by: T. Esposti Ongaro, C. Cavazzoni
!******************************************************************************
!
      PROGRAM pdac

      USE blunt_body, ONLY: set_blunt, ibl
      USE control_flags, ONLY: run
      USE dimensions
      USE domain_decomposition, ONLY: partition, ghost
      USE eos_gas, ONLY: allocate_eosg
      USE gas_constants, ONLY: allocate_gas_constants
      USE gas_solid_density, ONLY: allocate_density
      USE gas_solid_velocity, ONLY: allocate_velocity
      USE gas_solid_temperature, ONLY: allocate_temperature
      USE gas_solid_viscosity, ONLY: allocate_viscosity
      USE grid, ONLY: grid_setup, allocate_blbody, allocate_grid
      USE immersed_boundaries, ONLY: set_forcing, immb
      USE initial_conditions, ONLY: setpar, setup, cnvert, allocate_setup, npr
      USE input_module, ONLY: input, initc, number_of_block
      USE io_restart, ONLY: taperd, tapewr
      USE output_dump, ONLY: outp_recover
      USE parallel, ONLY: parallel_startup, parallel_hangup, &
     &    mpime, root, nproc
      USE particles_constants, ONLY: allocate_part_constants
      USE phases_matrix, ONLY: allocate_matrix
      USE pressure_epsilon, ONLY: allocate_press_eps
      USE specific_heat_module, ONLY: allocate_hcapgs
      USE tilde_momentum, ONLY: allocate_momentum
      USE time_parameters, ONLY: time, tstop, dt, tpr, tdump, itd, & 
     &                            timestart, rungekut
      USE turbulence_model, ONLY: allocate_turbo
      USE vent_conditions, ONLY: ivent, locate_vent
      USE dome_conditions, ONLY: idome, locate_dome
      USE volcano_topography, ONLY: import_topography, write_profile, itp
      USE environment, ONLY: cpclock, timing, elapsed_seconds
      USE io_files
!
      IMPLICIT NONE
      CHARACTER(LEN=3) :: procnum
      INTEGER :: mydate(10)
!
      INTEGER :: ig
      REAL*8 :: s0, s1, s2, s3, s4, s5, s6, t0
      REAL*8 :: pt0, pt1, pt2, pt3, pt4, pt5, pt6
      REAL*8 :: timtot, timprog, timdist, timsetup, timinit, timres, timghost
      REAL*8 :: mptimtot, mptimprog, mptimdist, mptimsetup, &
     &          mptiminit, mptimres, mptimghost
      LOGICAL :: debug = .FALSE.
      LOGICAL :: topen
!
! ... initialize parallel environment

      CALL parallel_startup  
!
! ... Initialize the system clock counter    
!
      t0 = elapsed_seconds()

          IF(timing) then
             s0 = cpclock()
             call MP_WALLTIME(pt0,mpime)
          END IF

! ... date and time
      CALL date_and_time( values = mydate )

      testnb = testfile//procnum(mpime)
      IF(mpime == root) THEN
        OPEN(UNIT=inputunit, FILE=inputfile, STATUS='UNKNOWN')
        IF(.NOT.debug ) OPEN(UNIT=logunit,FILE=logfile,STATUS='UNKNOWN')
        OPEN(UNIT=testunit,  FILE=testfile,  STATUS='UNKNOWN')
        OPEN(UNIT=errorunit, FILE=errorfile, STATUS='UNKNOWN')
        OPEN(UNIT=checkunit, FILE=checkfile, STATUS='UNKNOWN')
      ELSE
        OPEN(UNIT=testunit,  FILE=testnb,    STATUS='UNKNOWN')
      END IF

      IF( mpime == root ) THEN
        WRITE(logunit,100) mydate(5), mydate(6), mydate(7), mydate(3), mydate(2), mydate(1)
        WRITE(logunit,*)
        WRITE(logunit,*) 'Number of processor in use: ', nproc
      END IF
100   FORMAT( ' Pyroclastic Dispersion Analysis Code', /, &
           &  ' version: 3.0, June 2003',/, &
           &  ' authors: A.Neri, G.Macedonio, D.Gidaspow, T. Esposti Ongaro',/, &
           &  ' parallelized by: T. Esposti Ongaro, C. Cavazzoni, A. Neri',/, &
           &  ' 3d version implemented by: T. Esposti Ongaro, C. Cavazzoni',//, &
           &  ' This run began at ', I2, ':', I2, ':', I2, 3X, 'day ', I2, &
           &  ' month ', I2, ' year ', I4 )
!
! ... Read Input file
!
      CALL input( inputunit )
!
! ... set dimensions ...
!
      no = number_of_block
      nphase = nsolid + 1
!
! ... allocate global arrays 
!
      CALL allocate_grid
      CALL allocate_blbody
      CALL allocate_gas_constants
      CALL allocate_matrix
      CALL allocate_part_constants
      CALL allocate_setup
!
! ... initialize input fields
      CALL initc
!
! ... set start time
      timestart = time
!
! ... Setup cell-sizes and cell-flags
!
      CALL grid_setup

          IF (timing) then
              s1 = cpclock()
              call MP_WALLTIME(pt1,mpime)
          END IF   

! ... Import the topography from dem file and interpolate
! ... on the computational mesh
!
      IF (itp >= 1)  CALL import_topography

! ... Define volcanic vent position on the 3D topography
!
      IF (ivent >= 1) CALL locate_vent

! ... Define volcanic dome position on the 3D topography
!
      IF (idome >= 1) CALL locate_dome

! ... Set immersed boundary parameters (if prescribed)
!
      IF (immb == 1) CALL set_forcing

! ... Write the implicit topographic profile
!
      IF (itp >= 1)  CALL write_profile
!
      CALL partition

          IF (timing) then
              s2 = cpclock()
              call MP_WALLTIME(pt2,mpime)
          END IF
!
! ... Setting the ghost cells for parallel data exchange
! ... and the indexes
!
      CALL ghost

          IF (timing) then
              s3 = cpclock()
              call MP_WALLTIME(pt3,mpime)
          END IF
!
      IF (ibl >= 1) CALL set_blunt
!
! ... Attention!: all further operations are executed
! ... by each processor on local data
!
! ... Allocate local arrays
!
      CALL allocate_velocity
      CALL allocate_momentum
      CALL allocate_density
      CALL allocate_press_eps
      CALL allocate_temperature
      CALL allocate_eosg
      CALL allocate_viscosity
      CALL allocate_hcapgs
      CALL allocate_turbo
!
! ... Set parameters and initial conditions
!
      CALL setpar
      CALL setup

          IF (timing) then
              s4 = cpclock()
              call MP_WALLTIME(pt4,mpime)
          END IF
!
! ... Read restart file or recover initial
! ... conditions from an output file
!
      IF(itd == 2) THEN 
        CALL taperd
      ELSE IF (itd > 2) THEN
        CALL outp_recover
      END IF
!
! ... Compute initial conditions depending on restart mode
! ... Compute derived thermodynamic quantities.
!
      CALL cnvert

          IF (timing) then
              s5 = cpclock()
              call MP_WALLTIME(pt5,mpime)
          END IF
!
! ... Time advancement loop
!
STOP
      CALL prog
!
          IF (timing ) THEN
            s6 = cpclock()
            call MP_WALLTIME(pt6,mpime)
            timtot     = (s6 - s0)/1000.D0
            timprog    = (s6 - s5)/1000.D0
            timres     = (s5 - s4)/1000.D0
            timsetup   = (s4 - s3)/1000.D0
            timghost   = (s3 - s1)/1000.D0
            timinit    = (s1 - s0)/1000.D0
            mptimtot   = (pt6 - pt0)
            mptimprog  = (pt6 - pt5)          
            mptimres   = (pt5 - pt4)          
            mptimsetup = (pt4 - pt3)
            mptimghost = (pt3 - pt1)         
            mptiminit  = (pt1 - pt0)
             
          IF (mpime == root) THEN
            WRITE(logunit,*)' (From main) WALL TIME computed calling SYSTEM_CLOCK (s)'
            WRITE(logunit,900) 'Init', 'Ghost', 'Rest', 'Setup', 'Prog', 'Total'
            WRITE(logunit,999) timinit, timghost, timres, timsetup, timprog, timtot
!            WRITE(logunit,*)'             WALL TIME computed calling MP_WALLTIME (s)'
!            WRITE(logunit,900) 'Init', 'Ghost', 'Rest', 'Setup', 'Prog', 'Total'
!            WRITE(logunit,999) mptiminit, mptimghost, mptimres, mptimsetup, mptimprog, mptimtot
          END IF
999     FORMAT(6(1X,F10.2),/)
900     FORMAT(6(1X,A10))
          END IF

! ... date and time
      CALL date_and_time( values = mydate )
      IF( mpime == root ) THEN
        WRITE(logunit,110) mydate(5), mydate(6), mydate(7), mydate(3), mydate(2), mydate(1)
      END IF
110   FORMAT( ' This run ended at ', I2, ':', I2, ':', I2, 3X, 'day ', I2, &
              ' month ', I2, ' year ', I4 )
!
      INQUIRE(UNIT=inputunit,OPENED=topen)
      IF(topen) CLOSE(inputunit)
      IF( .NOT. debug ) THEN
        INQUIRE(UNIT=logunit,OPENED=topen)
        IF(topen) CLOSE(logunit)
      END IF
      INQUIRE(UNIT=errorunit,OPENED=topen)
      IF(topen) CLOSE(errorunit)
      INQUIRE(UNIT=testunit,OPENED=topen)
      IF(topen) CLOSE(testunit)
      INQUIRE(UNIT=checkunit,OPENED=topen)
      IF(topen) CLOSE(checkunit)
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
