!----------------------------------------------------------------------
      MODULE postp_input
!----------------------------------------------------------------------
      USE output_dump, ONLY: formatted_output
      IMPLICIT NONE
! ... read PP input file
! ... initialize some parameters
!
      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE postin( punit )

      USE filter_outp, ONLY: first_out, last_out, incr_out
      USE filter_outp, ONLY: downsize_x, downsize_y, downsize_z
      USE filter_outp, ONLY: act, number_of_points
      USE filter_outp, ONLY: assign_index
      USE filter_outp, ONLY: index_i, index_j, index_k
      USE filter_outp, ONLY: coordinate_x, coordinate_y, coordinate_z
      USE filter_outp, ONLY: variable_n, field_n
      USE process_outp, ONLY: imap, deltaz

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: punit
      
      NAMELIST / control / act, first_out, last_out, incr_out,  &
               downsize_x, downsize_y, downsize_z, formatted_output

      NAMELIST / map / imap, deltaz

      NAMELIST / sampling / number_of_points, assign_index, index_i, &
              index_j, index_k, coordinate_x, coordinate_y, coordinate_z

      NAMELIST / animation / variable_n

      NAMELIST / post_processing / field_n

!:::::::::::::::::::::::::::::  Sets default values  :::::::::::::::::::
!

! ... Control
      
  act        = 1        ! 1: frame processing, 2: time-space sampling        
  first_out  = 1        ! index of the first frame to be postprocessed
  last_out   = 10       ! index of the last frame to be postprocessed
  incr_out   = 1        ! increment between frame index
  downsize_x = 1    
  downsize_y = 1    
  downsize_z = 1    
  formatted_output = .TRUE.

! ... Map

  imap = 1
  deltaz = 20

! ... Sampling

  number_of_points = 1 
  assign_index     = .TRUE.
  index_i          = 1
  index_k          = 1
  coordinate_x     = 0.0
  coordinate_y     = 0.0
  coordinate_z     = 0.0

! ... Animation
  variable_n = 0

! ... Post Processing
  field_n = 0

!:::::::::::::::::::::::::::::  Read Namelists  ::::::::::::::::::::::::
!
      READ(punit, control)
      READ(punit, map)
      READ(punit, sampling)
      READ(punit, animation)
      READ(punit, post_processing)

      RETURN
      END SUBROUTINE postin
!----------------------------------------------------------------------
      END MODULE postp_input
!----------------------------------------------------------------------
