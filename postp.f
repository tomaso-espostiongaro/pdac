! ----------------------------------------------------------------------
      PROGRAM postp

      USE dimensions
      USE gas_constants, ONLY: allocate_gas_constants
      USE grid, ONLY: flic, allocate_blbody, allocate_grid, grid_setup, zzero
      USE initial_conditions, ONLY: allocate_setup, setpar
      USE input_module, ONLY: input, initc, number_of_block
      USE particles_constants, ONLY: allocate_part_constants
      USE postp_input, ONLY: postin
      USE process_outp, ONLY: filter, process
!
      IMPLICIT NONE
      CHARACTER(LEN=8) :: inputfile, logfile, errorfile, postfile
!
      INTEGER :: inputunit, logunit, errorunit, postunit
      INTEGER :: ig
!
! ... I/O files
!
      inputunit = 5
      logunit   = 6
      postunit  = 7
      errorunit = 8
      inputfile = 'pdac.dat'
      postfile = 'pp.dat'
      logfile = 'pp.log'
      errorfile = 'pp.err'

      OPEN(UNIT=inputunit, FILE=inputfile, STATUS='UNKNOWN')
      OPEN(UNIT=logunit,   FILE=logfile,   STATUS='UNKNOWN')
      OPEN(UNIT=postunit,  FILE=postfile,  STATUS='UNKNOWN')
      OPEN(UNIT=errorunit, FILE=errorfile, STATUS='UNKNOWN')

! ... Read Input files
!
      CALL input( inputunit )
      CALL postin( postunit )
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

      !CALL filter
      CALL process
!
      CLOSE(inputunit)
      CLOSE(logunit)
      CLOSE(postunit)
      CLOSE(errorunit)
!
      STOP
!
      END PROGRAM postp
!----------------------------------------------------------------------
