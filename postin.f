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

      USE filter_outp, ONLY: first_out, last_out, incr_out
      USE filter_outp, ONLY: downsize_x, downsize_y, downsize_z
      USE filter_outp, ONLY: act, number_of_probes, probe_file
      USE filter_outp, ONLY: assign_index
      USE filter_outp, ONLY: variable_n, field_n
      USE mass_orthoflux, ONLY: number_of_planes, planes_file
      USE mass_partition, ONLY: number_of_boxes, boxes_file
      USE process_outp, ONLY: imap, deltaz

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: punit
      
      NAMELIST / control / act, first_out, last_out, incr_out,  &
               downsize_x, downsize_y, downsize_z

      NAMELIST / map / imap, deltaz

      NAMELIST / sampling / number_of_probes, assign_index, probe_file

      NAMELIST / masspart / number_of_boxes, boxes_file

      NAMELIST / massflux / number_of_planes, planes_file

      NAMELIST / animation / variable_n

      NAMELIST / post_processing / field_n

!:::::::::::::::::::::::::::::  Sets default values  :::::::::::::::::::
!

! ... Control
      
  act        = 1        ! 0: filtering, 1: frame processing, 2: time-space sampling        
  first_out  = 0        ! index of the first frame to be postprocessed
  last_out   = 0        ! index of the last frame to be postprocessed
  incr_out   = 1        ! increment between frame index
  downsize_x = 1    
  downsize_y = 1    
  downsize_z = 1    

! ... Map

  imap = 1
  deltaz = 20

! ... Sampling

  number_of_probes = 1 
  probe_file       = 'probe.dat'
  assign_index     = .TRUE.

! ... Mass partition
  number_of_boxes = 1
  boxes_file = 'boxes.dat'

! ... Mass orthofluxes
  number_of_planes = 1
  planes_file = 'planes.dat'

! ... Animation
  variable_n = 0

! ... Post Processing
  field_n = 0

!:::::::::::::::::::::::::::::  Read Namelists  ::::::::::::::::::::::::
!
      READ(punit, control)
      READ(punit, map)
      READ(punit, sampling)
      READ(punit, masspart)
      READ(punit, massflux)
      READ(punit, animation)
      READ(punit, post_processing)

      RETURN
      END SUBROUTINE postin
!----------------------------------------------------------------------
      END MODULE postp_input
!----------------------------------------------------------------------
