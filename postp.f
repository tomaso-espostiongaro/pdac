!---------------------------------------------------------------------
! ... This program is a post-processing tool for PDAC outputs.
! ... It provides routines to filter large data-set and to compute
! ... derived physical fields from the primary fields computed by PDAC.
! ... It is intended for use on serial machines.
!
      PROGRAM postp

      USE dimensions
      USE filter_outp, ONLY: filter, read_implicit_profile, act
      USE gas_constants, ONLY: allocate_gas_constants
      USE grid, ONLY: flic, allocate_blbody, allocate_grid, grid_setup, zzero
      USE initial_conditions, ONLY: allocate_setup, setpar
      USE input_module, ONLY: input, initc, number_of_block
      USE io_files, ONLY: inputunit,postunit,logunit,testunit,ppunit,inputfile
      USE mass_orthoflux, ONLY: fluxn
      USE mass_partition, ONLY: massn
      USE parallel, ONLY: parallel_startup, parallel_hangup
      USE particles_constants, ONLY: allocate_part_constants
      USE postp_input, ONLY: postin
      USE postp_output, ONLY: write_avs_files, write_xml_files
      USE process_outp, ONLY: process, sample
      USE time_parameters, ONLY: itd
      USE volcano_topography, ONLY: itp
      USE control_flags, ONLY: prog
!
      IMPLICIT NONE
      CHARACTER(LEN=8) :: postfile, logfile, testfile, ppfile
!
! ... Define the parallel environment
! ... (needed to use routines written for PDAC)
!

      prog = 'POSP'

      CALL parallel_startup
!
! ... Default I/O files
!
      logfile   = 'pp.log'
      testfile  = 'pp.tst'
      postfile  = 'pp.dat'
      ppfile    = 'dpd.dat'
!
      OPEN(UNIT=inputunit, FILE=inputfile, STATUS='UNKNOWN')
      OPEN(UNIT=postunit,  FILE=postfile,  STATUS='UNKNOWN')
      OPEN(UNIT=logunit, FILE=logfile, STATUS='UNKNOWN')
      OPEN(UNIT=testunit,  FILE=testfile,  STATUS='UNKNOWN')
      OPEN(UNIT=ppunit,  FILE=ppfile,  STATUS='UNKNOWN')
!
! ... Read Input files 'pdac.dat'
!
      CALL input( inputunit )
!
! ... Read Input files
!
      CALL postin( postunit )
!
! ... By default, postp cannot restart
      itd = 1
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
!
! ... Initialize some input variables
!
      CALL initc
!
! ... Set cell-type flags
!
      CALL grid_setup
!
! ... Read the implicit profile
!
      IF (itp >= 1) CALL read_implicit_profile
!
! ... Set physical parameters and useful constants
!
      CALL setpar
!
! ... Here start the post-processing core
!********************************************************
!
! ... Split OUTPUT files, downsize, crop, etc.
!
      IF (act == 0) CALL filter
!
! ... Write input files for visualization tools
!
      CALL write_avs_files
      CALL write_xml_files
!
! ... Compute derived fields from the primary OUTPUT fields;
! ... map hazard variables at a given height above ground
!
      IF (act == 1) CALL process
      IF (act == 2) CALL sample
      IF (act == 3) CALL massn
      IF (act == 4) CALL fluxn
!
!********************************************************
!
      CLOSE(inputunit)
      CLOSE(postunit)
      CLOSE(logunit)
      CLOSE(testunit)
      CLOSE(ppunit)

! ... close parallel environment
!
      CALL parallel_hangup

      STOP
      END PROGRAM postp
!----------------------------------------------------------------------
