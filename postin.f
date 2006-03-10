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

      USE filter_outp, ONLY: downsize_x, downsize_y, downsize_z
      USE filter_outp, ONLY: variable_n, field_n
      USE sample_points, ONLY: number_of_probes, assign_index, probe_file, isamp
      USE mass_orthoflux, ONLY: number_of_planes, planes_file, ifluxn
      USE mass_partition, ONLY: number_of_boxes, boxes_file, imassn
      USE mass_ground, ONLY: thickness, iground
      USE process_outp, ONLY: act
      USE process_outp, ONLY: iflds, imap
      USE io_serial, ONLY: first_out, last_out, incr_out
      USE io_serial, ONLY: deltaz

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: punit
      
      NAMELIST / control / act, first_out, last_out, incr_out,  &
               downsize_x, downsize_y, downsize_z

      NAMELIST / fields / iflds

      NAMELIST / map / imap, deltaz

      NAMELIST / sampling / isamp, number_of_probes, assign_index, probe_file

      NAMELIST / masspart / imassn, number_of_boxes, boxes_file

      NAMELIST / massflux / ifluxn, number_of_planes, planes_file

      NAMELIST / massgsedim / iground, thickness

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

! ... Fields
  
  iflds = 0

! ... Map

  imap = 1
  deltaz = 20

! ... Sampling
  
  isamp = 0
  number_of_probes = 1 
  probe_file       = 'probe.dat'
  assign_index     = .TRUE.

! ... Mass partition
 
  imassn = 0
  number_of_boxes = 1
  boxes_file = 'boxes.dat'

! ... Mass orthofluxes

  ifluxn = 0
  number_of_planes = 1
  planes_file = 'planes.dat'

! ... Mass ground deposit
  
  iground = 0
  thickness = 500.D0

! ... Animation
  variable_n = 0

! ... Post Processing
  field_n = 0

!:::::::::::::::::::::::::::::  Read Namelists  ::::::::::::::::::::::::
!
      READ(punit, control)
      READ(punit, fields)
      READ(punit, map)
      READ(punit, sampling)
      READ(punit, masspart)
      READ(punit, massflux)
      READ(punit, massgsedim)
      READ(punit, animation)
      READ(punit, post_processing)

      RETURN
      END SUBROUTINE postin
!----------------------------------------------------------------------
      END MODULE postp_input
!----------------------------------------------------------------------
