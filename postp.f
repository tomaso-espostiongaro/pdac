!----------------------------------------------------------------------
      PROGRAM postp

      USE dimensions
      USE gas_constants, ONLY: allocate_gas_constants
      USE grid, ONLY: flic, allocate_blbody, allocate_grid, grid_setup, zzero
      USE initial_conditions, ONLY: allocate_setup, setpar
      USE input_module, ONLY: input, initc, number_of_block
      USE particles_constants, ONLY: allocate_part_constants
      USE postp_input, ONLY: postin
      USE postp_output, ONLY: write_avs_files, write_xml_files
      USE process_outp, ONLY: filter, process
      USE parallel, ONLY: parallel_startup, parallel_hangup, &
     &    mpime, root, nproc
      USE io_files, ONLY: inputunit, postunit, logunit, testunit, ppunit, inputfile

!
      IMPLICIT NONE
      CHARACTER(LEN=8) :: inputfile, postfile, logfile, testfile, ppfile
!
      INTEGER :: ig

!
! ... Initialize parallel environment
!

      CALL parallel_startup

!
! ... I/O files
!

      logfile   = 'pp.log'
      testfile  = 'pp.tst'
      postfile  = 'pp.dat'
      ppfile    = 'dpd.dat'

      IF( mpime == 0 ) THEN

        OPEN(UNIT=inputunit, FILE=inputfile, STATUS='UNKNOWN')
        OPEN(UNIT=postunit,  FILE=postfile,  STATUS='UNKNOWN')
        OPEN(UNIT=logunit, FILE=logfile, STATUS='UNKNOWN')
        OPEN(UNIT=testunit,  FILE=testfile,  STATUS='UNKNOWN')
        OPEN(UNIT=ppunit,  FILE=ppfile,  STATUS='UNKNOWN')

      END IF

!
! ... Read Input files
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
      CALL allocate_part_constants
      CALL allocate_setup

      CALL initc
!
! ... Set cell-type flags
!
      CALL flic
!
      CALL grid_setup
!
      CALL setpar

! ... Read Input files
!
      CALL postin( postunit )

      !CALL filter

      CALL write_avs_files
      CALL write_xml_files

      CALL process
!
      IF( mpime == 0 ) THEN
        CLOSE(inputunit)
        CLOSE(postunit)
        CLOSE(logunit)
        CLOSE(testunit)
        CLOSE(ppunit)
      END IF

      CALL parallel_hangup
!
      STOP
!
      END PROGRAM postp
!----------------------------------------------------------------------
