!----------------------------------------------------------------------
      PROGRAM postp

      USE dimensions
      USE gas_constants, ONLY: allocate_gas_constants
      USE grid, ONLY: flic, allocate_blbody, allocate_grid, grid_setup, zzero
      USE initial_conditions, ONLY: allocate_setup, setpar
      USE input_module, ONLY: input, initc, number_of_block
      USE particles_constants, ONLY: allocate_part_constants
      USE postp_input, ONLY: postin
      USE postp_output, ONLY: write_avs_files
      USE process_outp, ONLY: filter, process
!
      IMPLICIT NONE
      CHARACTER(LEN=8) :: inputfile, postfile
!
      INTEGER :: inputunit, postunit
      INTEGER :: ig
!
! ... I/O files
!
      inputunit = 5
      postunit  = 7
      inputfile = 'pdac.dat'
      postfile = 'pp.dat'

      OPEN(UNIT=inputunit, FILE=inputfile, STATUS='UNKNOWN')
      OPEN(UNIT=postunit,  FILE=postfile,  STATUS='UNKNOWN')

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

      CALL write_avs_files

! ... Read Input files
!
      CALL postin( postunit )
      !CALL filter
      CALL process
!
      CLOSE(inputunit)
      CLOSE(postunit)
!
      STOP
!
      END PROGRAM postp
!----------------------------------------------------------------------
