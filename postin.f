!----------------------------------------------------------------------
      MODULE postp_input
!----------------------------------------------------------------------
      IMPLICIT NONE
! ... read PP input file
! ... initialize some parameters
!
      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE postin( punit )

      USE process_outp, ONLY: first_out, last_out, incr_out
      USE process_outp, ONLY: downsize_x, downsize_y, downsize_z
      USE process_outp, ONLY: act, number_of_points
      USE process_outp, ONLY: assign_index
      USE process_outp, ONLY: index_i, index_j, index_k
      USE process_outp, ONLY: coordinate_x, coordinate_y, coordinate_z
      USE process_outp, ONLY: variable_n, field_n
      USE process_outp, ONLY: formatted_output

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: punit
      
      NAMELIST / control / act, first_out, last_out, incr_out,  &
               downsize_x, downsize_y, downsize_z, formatted_output

      NAMELIST / sampling / number_of_points, assign_index, index_i, &
              index_j, index_k, coordinate_x, coordinate_y, coordinate_z

      NAMELIST / animation / variable_n

      NAMELIST / post_processing / field_n

!:::::::::::::::::::::::::::::  Sets default values  :::::::::::::::::::
!

! ... Control
      
  act        = 2        ! 1: mount time sequence, 2: time-space sampling        
  first_out  = 1     
  last_out   = 10   
  incr_out   = 1    
  downsize_x = 1    
  downsize_y = 1    
  downsize_z = 1    
  formatted_output = .TRUE.

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
      READ(punit, sampling)
      READ(punit, animation)
      READ(punit, post_processing)

      RETURN
      END SUBROUTINE postin
!----------------------------------------------------------------------
      END MODULE postp_input
!----------------------------------------------------------------------
