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

#if defined __HPM
#  include "/cineca/prod/hpm/include/f_hpm.h"
#endif
!
!
      PROGRAM pp

      USE dimensions
      USE domain_decomposition, ONLY: partition, ghost
      USE eos_gas, ONLY: allocate_eosg
      USE gas_constants, ONLY: allocate_gas_constants
      USE gas_solid_density, ONLY: allocate_density
      USE gas_solid_velocity, ONLY: allocate_velocity
      USE gas_solid_temperature, ONLY: allocate_temperature
      USE gas_solid_viscosity, ONLY: allocate_viscosity
      USE grid, ONLY: flic, allocate_blbody, allocate_grid, grid_setup
      USE initial_conditions, ONLY: allocate_setup, zzero, setc
      USE input_module, ONLY: input, initc, number_of_block
      USE input_module, ONLY: first_out, last_out, incr_out
      USE parallel, ONLY: parallel_startup, parallel_hangup, &
     &    mpime, root
      USE particles_constants, ONLY: allocate_part_constants
      USE phases_matrix, ONLY: allocate_matrix
      USE pressure_epsilon, ONLY: allocate_press_eps
      USE roughness_module, ONLY: zrough, deallocate_roughness
      USE specific_heat_module, ONLY: allocate_hcapgs
      USE tilde_momentum, ONLY: allocate_momentum
      USE time_parameters, ONLY: time, tstop, dt, tpr, tdump, itd, & 
     &                            timestart, rungekut
      USE turbulence_model, ONLY: allocate_turbo
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
      logfile = 'pp.log'
      testfile = 'pp.tst'
      errorfile = 'pp.err'
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
      CALL input( inputunit, 'PP' )
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
      CALL allocate_part_constants
      CALL allocate_setup

      CALL initc
!
! ... Set cell-type flags
!
      CALL flic
!
! ... Domain decomposition for parallelization 
!
      CALL partition

!
! ... Setting ghost cells for parallel data exchange
! ... and the indexes
!
      CALL ghost
!
      CALL grid_setup( zzero )
      !CALL setc

      CALL filter( first_out, last_out, incr_out )
!
! ... terminate the IBM HW performance monitor session
      !call f_hpmterminate( mpime )

      CLOSE(5)
      IF( .NOT. debug ) CLOSE(6)
      CLOSE(7)
      CLOSE(8)
!
! ... Finalize parallel environment

      CALL parallel_hangup
!
      STOP
!
!**************************************************************
      END PROGRAM pp
!**************************************************************
