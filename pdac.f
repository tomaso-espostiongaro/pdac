!******************************************************************************
!> \mainpage
!> Parallel multiphase and multicomponent fluid-dynamic code 
!> for the simulation of pyroclastic dispersion processes
!>
!> \version 4.0
!> \date: Jan 2014
!> \author T. Esposti Ongaro, A.Neri, G.Macedonio, D.Gidaspow
!<******************************************************************************
   PROGRAM pdac
!
      USE blunt_body, ONLY: set_blunt, ibl
      USE compute_mean_fields, ONLY: allocate_mean_fields, imrt
      USE control_flags, ONLY: prog, lpr, nfil
      USE dimensions
      USE domain_decomposition, ONLY: partition
      USE domain_mapping, ONLY: ghost
      USE dome_conditions, ONLY: idome, locate_dome
      USE mixture_fields, ONLY: allocate_mixture_fields
      USE environment, ONLY: cpclock, timing, elapsed_seconds
      USE eos_gas, ONLY: allocate_eosg
      USE gas_components, ONLY: allocate_species
      USE gas_constants, ONLY: allocate_gas_constants
      USE gas_solid_density, ONLY: allocate_density
      USE gas_solid_velocity, ONLY: allocate_velocity
      USE gas_solid_temperature, ONLY: allocate_temperature
      USE gas_solid_viscosity, ONLY: allocate_viscosity
      USE grid, ONLY: grid_setup, allocate_blbody, allocate_grid
      USE immersed_boundaries, ONLY: set_forcing, immb
      USE initial_conditions, ONLY: setpar, setup, cnvert, allocate_setup
      USE input_module, ONLY: input, initc
      USE io_files, ONLY: logunit, dataunit, errorunit, checkunit, inputunit, testunit
      USE io_files, ONLY: logfile, datafile, errorfile, checkfile, inputfile, testfile
      USE io_files, ONLY: testnb
      USE io_restart, ONLY: taperd, tapewr, outp_recover, outp_remap
      USE mass_sink, ONLY: allocate_sink, isink, read_mass_loss
      USE parallel, ONLY: parallel_startup, parallel_hangup, &
     &    mpime, root, nproc
      USE particles_constants, ONLY: allocate_part_constants
      USE phases_matrix, ONLY: allocate_matrix
      USE pressure_epsilon, ONLY: allocate_press_eps
      USE specific_heat_module, ONLY: allocate_hcapgs
      USE tilde_momentum, ONLY: allocate_momentum
      USE time_advancement, ONLY: main_prog
      USE time_parameters, ONLY: time, tstop, dt, tpr, tdump, itd, & 
     &                            timestart, rungekut
      USE turbulence_model, ONLY: allocate_turbo
      USE vent_conditions, ONLY: ivent, locate_vent
      USE volcano_topography, ONLY: import_topography, write_profile, itp, &
                                    export_topography
!
      IMPLICIT NONE
      CHARACTER(LEN=3) :: procnum
      INTEGER :: mydate(10)
!
      INTEGER :: ig, n
      REAL*8 :: s0, s1, s2, s3, s4, s5, s6, t0, tt
      REAL*8 :: pt0, pt1, pt2, pt3, pt4, pt5, pt6
      REAL*8 :: timtot, timprog, timdist, timsetup, timinit, timres, timghost
      REAL*8 :: mptimtot, mptimprog, mptimdist, mptimsetup, &
     &          mptiminit, mptimres, mptimghost
      LOGICAL :: topen
      INTEGER :: ios
!
! ... initialize parallel environment

      prog = 'PDAC'

      CALL parallel_startup  
!
! ... Initialize the system clock counter    
!
      t0 = elapsed_seconds()

          IF(timing) then
             s0 = cpclock()
             !call MP_WALLTIME(pt0,mpime)
          END IF
!
! ... date and time
!
      CALL date_and_time( values = mydate )
!
! ... Root processor opens log and error files
!
      IF(mpime == root) THEN
        OPEN(UNIT=logunit,FILE=logfile,STATUS='UNKNOWN')
        OPEN(UNIT=dataunit,FILE=datafile,STATUS='UNKNOWN')
        OPEN(UNIT=errorunit, FILE=errorfile, STATUS='UNKNOWN')
        OPEN(UNIT=checkunit, FILE=checkfile, POSITION='APPEND',STATUS='UNKNOWN')
      END IF
!
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
      IF(mpime == root) OPEN(UNIT=inputunit, FILE=inputfile, STATUS='UNKNOWN',ERR=199)
      CALL input( inputunit )
      IF(mpime == root) CLOSE(inputunit)
!
! ... Open Test files 
! ... (WARNING!: opening 'nproc' files on the same unit could cause 
! ... file-system problems on some architectures)
!
      testnb = testfile//procnum(mpime)
      IF (lpr > 0) THEN
        IF(mpime == root) THEN
          OPEN(UNIT=testunit,  FILE=testfile,  STATUS='UNKNOWN')
        ELSE
          OPEN(UNIT=testunit,  FILE=testnb,    STATUS='UNKNOWN')
        END IF
      END IF
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
! ... Setup cell-sizes and cell flags
!
      CALL grid_setup
!
          IF (timing) then
              s1 = cpclock()
              !call MP_WALLTIME(pt1,mpime)
          END IF   

! ... Import the topography from dem file and interpolate
! ... on the computational mesh. Change the cell flags.
!
      IF (itp >= 1)  CALL import_topography

          IF (timing) then
              s2 = cpclock()
              IF( mpime == root ) WRITE(logunit,*) 'time for import_topography = ', s2-s1
              tt = s2
          END IF

! ... Define volcanic vent position on the 3D topography
! ... and change the cell flags (can modify the topography)
!
      IF (ivent >= 1) CALL locate_vent

          IF (timing) then
              s2 = cpclock()
              IF( mpime == root ) WRITE(logunit,*) 'time for locate_vent = ', s2-tt
              tt = s2
          END IF

! ... Define volcanic dome position on the 3D topography
! ... and change the cell flags (can modify the topography)
!
      IF (idome >= 1) CALL locate_dome

          IF (timing) then
              s2 = cpclock()
              IF( mpime == root ) WRITE(logunit,*) 'time for locate_dome = ', s2-tt
              tt = s2
          END IF

! ... Initialize the array of elevations on the computational mesh
! 
      CALL export_topography

          IF (timing) then
              s2 = cpclock()
              IF( mpime == root ) WRITE(logunit,*) 'time for export_topography = ', s2-tt
              tt = s2
          END IF

! ... Set immersed boundary parameters (if prescribed)
! ... and change the cell flags
!
      IF (immb == 1) CALL set_forcing

          IF (timing) then
              s2 = cpclock()
              IF( mpime == root ) WRITE(logunit,*) 'time for set_forcing = ', s2-tt
              tt = s2
          END IF

! ... Write the implicit topographic profile
!
      IF (itp >= 1) CALL write_profile

          IF (timing) then
              s2 = cpclock()
              IF( mpime == root ) WRITE(logunit,*) 'time for write_profile = ', s2-tt
              tt = s2
          END IF
!
      IF (itd == -1) THEN
        CALL parallel_hangup
        STOP
      END IF
!
! ... Domain Decomposition
!
      !Simulate the domain decomposition!
      !DO n = 512, 512
      !CALL partition( n, mpime, root ) ! TEST
      !END DO
      !CALL parallel_hangup   ! TEST
      !STOP                   ! TEST
      !
      CALL partition( nproc, mpime, root )

          IF (timing) then
              s2 = cpclock()
              !call MP_WALLTIME(pt2,mpime)
              IF( mpime == root ) WRITE(logunit,*) 'time for partition = ', s2-tt
          END IF
!
! ... Setting the ghost cells for parallel data exchange
! ... and the indexes. Change the cell flags.
!
      CALL ghost

          IF (timing) then
              s3 = cpclock()
              !call MP_WALLTIME(pt3,mpime)
              IF( mpime == root ) WRITE(logunit,*) 'time for ghost = ', s3-s2
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
      CALL allocate_species
      CALL allocate_press_eps
      CALL allocate_temperature
      CALL allocate_eosg
      CALL allocate_viscosity
      CALL allocate_hcapgs
      CALL allocate_turbo
      CALL allocate_sink
      CALL allocate_mixture_fields
      IF (imrt >= 1) CALL allocate_mean_fields
!
! ... Set parameters and initial conditions
!
      CALL setpar

      CALL setup

          IF (timing) then
              s4 = cpclock()
              !call MP_WALLTIME(pt4,mpime)
          END IF
!
! ... Read restart file or recover initial
! ... conditions from an output file
!
      IF(itd == 2) THEN 
        IF (isink > 0) CALL read_mass_loss
        CALL taperd
      ELSE IF (itd == 3) THEN
        IF (isink > 0) CALL read_mass_loss(nfil)
        CALL outp_recover(nfil)
      ELSE IF (itd == 4) THEN
        CALL outp_remap(nfil)
      END IF
!
! ... reset start time
      timestart = time
!
! ... Set the initial conditions in the ghost cells.
! ... Compute initial conditions depending on restart mode.
! ... Compute derived thermodynamic quantities.
!
      CALL cnvert
!
          IF (timing) then
              s5 = cpclock()
              !call MP_WALLTIME(pt5,mpime)
          END IF
!
! ... Time advancement loop
!
      CALL main_prog
!
          IF (timing ) THEN
            s6 = cpclock()
            !call MP_WALLTIME(pt6,mpime)
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
      INQUIRE(UNIT=logunit,OPENED=topen)
      IF(topen) CLOSE(logunit)
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
 199  CALL error('main','error opening inputunit',inputunit)      
!
!**************************************************************
   END PROGRAM pdac
!**************************************************************
