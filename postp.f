!---------------------------------------------------------------------
! ... This program is a post-processing tool for PDAC outputs.
! ... It provides routines to filter large data-set and to compute
! ... derived physical fields from the primary fields computed by PDAC.
! ... It is intended for use on serial machines.
!
      PROGRAM postp

      USE control_flags, ONLY: prog, lpr
      USE dimensions
      USE domain_decomposition, ONLY: partition
      USE domain_mapping, ONLY: ghost
      USE filter_outp, ONLY: filter
      USE gas_constants, ONLY: allocate_gas_constants
      USE grid, ONLY: flic, allocate_blbody, allocate_grid, grid_setup, zzero
      USE initial_conditions, ONLY: allocate_setup, setpar
      USE input_module, ONLY: input, initc, number_of_block
      USE io_files, ONLY: inputunit,postunit,logunit,testunit, &
                          inputfile, errorunit
      USE parallel, ONLY: parallel_startup, parallel_hangup, mpime, root, nproc
      USE particles_constants, ONLY: allocate_part_constants
      USE postp_input, ONLY: postin
      USE postp_output, ONLY: read_implicit_profile
      USE process_outp, ONLY: process, act
      USE time_parameters, ONLY: itd
      USE volcano_topography, ONLY: itp
!
      IMPLICIT NONE
      CHARACTER(LEN=3) :: procnum
      CHARACTER(LEN=6) :: postfile, logfile, testfile, errfile
      CHARACTER(LEN=9) :: testnb
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
      errfile   = 'pp.err'
!
      IF (mpime == root) THEN
        OPEN(UNIT=inputunit, FILE=inputfile, STATUS='UNKNOWN')
        OPEN(UNIT=postunit,  FILE=postfile,  STATUS='UNKNOWN')
        OPEN(UNIT=logunit, FILE=logfile, STATUS='UNKNOWN')
        OPEN(UNIT=errorunit,  FILE=errfile,  STATUS='UNKNOWN')
        WRITE(logunit,*) 'Number of processor in use: ', nproc
      END IF
!
! ... Processor 'root' read Input file 'pdac.dat'
! ... and broadcasts data
!
      CALL input( inputunit )
      IF (mpime == root) WRITE(logunit,*) 'input file 1 read'
!
! ... Open Test files
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
! ... Processor 'root' read Input file 'pp.dat'
! ... and broadcasts data
!
      CALL postin( postunit )
      IF (mpime == root) WRITE(logunit,*) 'input file 2 read'
!
! ... By default, postp cannot restart
      itd = 1
!
! ... set dimensions ...
!
      no = number_of_block
      nphase = nsolid + 1
!
! ... allocate global arrays (all processes)
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
      IF (act == 0) GOTO 111 
!
! ... Set physical parameters and useful constants (All processes)
!
      CALL setpar
!
! ... Domain Decomposition (WARNING: the number of cells
! ... within subdomain could be unbalanced, due to the
! ... lack of the topography setup)
!
      CALL partition( nproc, mpime, root )
!
! ... Setting the ghost cells for parallel data exchange
! ... and the indexes
!
      CALL ghost
!
! ... Here start the post-processing core
!********************************************************
!
! ... Split OUTPUT files, downsize, crop, etc.
!
      ! CALL filter
!
! ... Compute derived fields from the primary OUTPUT fields;
! ... map hazard variables at a given height above ground
!
      IF (act == 1) CALL process
!
!********************************************************
!
 111  CONTINUE
      IF (mpime == root) THEN
        CLOSE(inputunit)
        CLOSE(postunit)
        CLOSE(logunit)
        CLOSE(testunit)
      END IF
!
! ... close parallel environment
!
      CALL parallel_hangup
      STOP
      END PROGRAM postp
!----------------------------------------------------------------------
