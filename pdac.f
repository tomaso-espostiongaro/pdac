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
!#if defined __HPM
!#  include "/cineca/prod/hpm/include/f_hpm.h"
!#endif
      PROGRAM pdac

      USE control_flags, ONLY: itp
      USE blunt_body, ONLY: set_blunt, ibl
      USE dimensions
      USE domain_decomposition, ONLY: partition, ghost
      USE eos_gas, ONLY: allocate_eosg
      USE gas_constants, ONLY: allocate_gas_constants
      USE gas_solid_density, ONLY: allocate_density
      USE gas_solid_velocity, ONLY: allocate_velocity
      USE gas_solid_temperature, ONLY: allocate_temperature
      USE gas_solid_viscosity, ONLY: allocate_viscosity
      USE grid, ONLY: grid_setup, allocate_blbody, allocate_grid
      USE immersed_boundaries, ONLY: import_topo
      USE initial_conditions, ONLY: setpar, setup, resetup, allocate_setup, npr
      USE input_module, ONLY: input, initc, number_of_block
      USE io_restart, ONLY: taperd, tapewr
      USE parallel, ONLY: parallel_startup, parallel_hangup, &
     &    mpime, root
      USE particles_constants, ONLY: allocate_part_constants
      USE phases_matrix, ONLY: allocate_matrix
      USE pressure_epsilon, ONLY: allocate_press_eps
      USE specific_heat_module, ONLY: allocate_hcapgs
      USE tilde_momentum, ONLY: allocate_momentum
      USE time_parameters, ONLY: time, tstop, dt, tpr, tdump, itd, & 
     &                            timestart, rungekut
      USE turbulence_model, ONLY: allocate_turbo
      USE vent_conditions, ONLY: ivent, locate_vent
      USE volcano_topography, ONLY: read_topo
      USE environment, ONLY: cpclock, timing, elapsed_seconds
!
      IMPLICIT NONE
      CHARACTER(LEN=11) :: errnb, testnb, lognb
      CHARACTER(LEN=8) :: inputfile, logfile, errorfile, testfile
      CHARACTER(LEN=3) :: procnum
      INTEGER :: mydate(10)
!
      INTEGER :: inputunit, logunit, errorunit, testunit
      INTEGER :: ig
      REAL*8 :: s0, s1, s2, s3, s4, s5, s6, t0
      REAL*8 :: pt0, pt1, pt2, pt3, pt4, pt5, pt6
      REAL*8 :: timtot, timprog, timdist, timsetup, timinit, timres, timghost
      REAL*8 :: mptimtot, mptimprog, mptimdist, mptimsetup, &
     &          mptiminit, mptimres, mptimghost
      LOGICAL :: debug = .FALSE.
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

! ... Initialize the IBM HW performance monitor
      !call f_hpminit( mpime, 'pdac' )
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
      lognb = logfile//procnum(mpime)
      IF(mpime .EQ. root) THEN
        OPEN(UNIT=inputunit, FILE=inputfile, STATUS='UNKNOWN')
        IF( .NOT. debug ) OPEN(UNIT=logunit,   FILE=logfile,   STATUS='UNKNOWN')
        OPEN(UNIT=testunit,  FILE=testfile,  STATUS='UNKNOWN')
        OPEN(UNIT=errorunit, FILE=errorfile, STATUS='UNKNOWN')
      ELSE
        IF( .NOT. debug ) OPEN(UNIT=logunit,   FILE=lognb,   STATUS='UNKNOWN')
        OPEN(UNIT=testunit,  FILE=testnb,    STATUS='UNKNOWN')
        OPEN(UNIT=errorunit, FILE=errnb,     STATUS='UNKNOWN')
      END IF

      WRITE(6,100) mydate(5), mydate(6), mydate(7), mydate(3), mydate(2), mydate(1)
100   FORMAT( ' Pyroclastic Dispersion Analysis Code', /, &
              ' version: 3.0, June 2003',/, &
              ' authors: A.Neri, G.Macedonio, D.Gidaspow, T. Esposti Ongaro',/, &
              ' parallelized by: T. Esposti Ongaro, C. Cavazzoni, A. Neri',/, &
              ' 3d version implemented by: T. Esposti Ongaro, C. Cavazzoni',//, &
              ' This run began at ', I2, ':', I2, ':', I2, 3X, 'day ', I2, &
              ' month ', I2, ' year ', I4 )
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

! ... Read topography file 
!
      IF (itp >= 1)  CALL read_topo

! ... Import the topography on the computational mesh 
! ... and set immersed boundary parameters
!
      IF (itp >= 1) CALL import_topo
 
! ... Define volcanic vent position on the 3D topography
!
      IF (ivent >= 1) CALL locate_vent
!
! ... Domain decomposition for parallelization 
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
! ... Read restart file
!
      IF(itd == 2) THEN 
        CALL taperd
      ELSE IF (itd > 2) THEN
        CALL error('setup','Output recovering not implemented',1)         
      END IF
!
! ... Re-compute initial conditions depending on restart mode
! ... (i.e. when itd > 2 )
!
      CALL resetup

      IF (timing) then
          s5 = cpclock()
          call MP_WALLTIME(pt5,mpime)
      END IF
!
! ... Time advancement loop
!
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
         
          WRITE(7,*)' (From main) WALL TIME computed calling SYSTEM_CLOCK (s)'
          WRITE(7,900) 'Init', 'Ghost', 'Rest', 'Setup', 'Prog', 'Total'
          WRITE(7,999) timinit, timghost, timres, timsetup, timprog, timtot
          WRITE(7,*)'             WALL TIME computed calling MP_WALLTIME (s)'
          WRITE(7,900) 'Init', 'Ghost', 'Rest', 'Setup', 'Prog', 'Total'
          WRITE(7,999) mptiminit, mptimghost, mptimres, mptimsetup, mptimprog, mptimtot
999       FORMAT(6(1X,F10.2),/)
900       FORMAT(6(1X,A10))
        END IF

! ... terminate the IBM HW performance monitor session
      !call f_hpmterminate( mpime )

! ... date and time
      CALL date_and_time( values = mydate )
      WRITE(6,110) mydate(5), mydate(6), mydate(7), mydate(3), mydate(2), mydate(1)
110   FORMAT( ' This run ended at ', I2, ':', I2, ':', I2, 3X, 'day ', I2, &
              ' month ', I2, ' year ', I4 )
!
      CLOSE(inputunit)
      IF( .NOT. debug ) CLOSE(logunit)
      CLOSE(errorunit)
      CLOSE(testunit)
      CLOSE(15)
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
