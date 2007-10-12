!----------------------------------------------------------------------
      MODULE postp_input
!----------------------------------------------------------------------
      USE parallel, ONLY: mpime, root
      IMPLICIT NONE
! ... read PP input file
! ... initialize some parameters
!
      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE postin( punit )

      USE dimensions, ONLY: max_size
      USE filter_outp, ONLY: downsize_x, downsize_y, downsize_z
      USE filter_outp, ONLY: variable_n, field_n
      USE sample_points, ONLY: number_of_probes, assign_index, probe_file, isamp
      USE mass_orthoflux, ONLY: number_of_planes, planes_file, ifluxn
      USE mass_partition, ONLY: number_of_boxes, boxes_file, imassn
      USE mass_ground, ONLY: thickness, iground
      USE process_outp, ONLY: act
      USE process_outp, ONLY: iflds, imnfld, imap
      USE postp_output, ONLY: first_out, last_out, incr_out
      USE postp_output, ONLY: iminc, imaxc, jminc, jmaxc, kminc, kmaxc 
      USE postp_output, ONLY: deltaz
      USE postp_output, ONLY: print_log, print_tg, print_mn, print_cm, print_pd, print_mnn 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: punit
      
      NAMELIST / control / act, first_out, last_out, incr_out,  &
               downsize_x, downsize_y, downsize_z, iminc, imaxc, jminc, &
               jmaxc, kminc, kmaxc

      NAMELIST / fields / iflds, print_log, print_tg, print_mn, print_cm, print_pd, print_mnn 

      NAMELIST / mean_outp / imnfld

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
  first_out  = 601        ! index of the first frame to be postprocessed
  last_out   = 710        ! index of the last frame to be postprocessed
  incr_out   = 1        ! increment between frame index
  downsize_x = 1    
  downsize_y = 1    
  downsize_z = 1    
  iminc = 1
  imaxc = max_size
  jminc = 1
  jmaxc = max_size
  kminc = 1
  kmaxc = max_size

! ... Fields
  
  iflds = 1
  print_log = .FALSE.
  print_tg  = .TRUE.
  print_mn  = .FALSE.
  print_cm  = .FALSE.
  print_pd  = .FALSE.
  print_mnn = .FALSE.

! ... Mean Field

  imnfld = 0

! ... Map

  imap = 0
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
! ... Control namelist .................................................
!
      IF (mpime == root) READ(punit, control)
      CALL bcast_integer(act,1,root)
      CALL bcast_integer(first_out,1,root)
      CALL bcast_integer(last_out,1,root)
      CALL bcast_integer(incr_out,1,root)
      CALL bcast_integer(downsize_x,1,root)
      CALL bcast_integer(downsize_y,1,root)
      CALL bcast_integer(downsize_z,1,root)
      CALL bcast_integer(iminc,1,root)
      CALL bcast_integer(imaxc,1,root)
      CALL bcast_integer(jminc,1,root)
      CALL bcast_integer(jmaxc,1,root)
      CALL bcast_integer(kminc,1,root)
      CALL bcast_integer(kmaxc,1,root)
!
! ... Fields namelist ..................................................
!
      IF (mpime == root) READ(punit, fields)
      CALL bcast_integer(iflds,1,root)
      CALL bcast_logical(print_log,1,root)
      CALL bcast_logical(print_tg,1,root)
      CALL bcast_logical(print_mn,1,root)
      CALL bcast_logical(print_cm,1,root)
      CALL bcast_logical(print_pd,1,root)
      CALL bcast_logical(print_mnn,1,root)
!
! ... Maen Field namelist ...............................................
!
      IF (mpime == root) READ(punit, mean_outp)
      CALL bcast_integer(imnfld,1,root)
!
! ... Map namelist ......................................................
!
      IF (mpime == root) READ(punit, map)
      CALL bcast_integer(imap,1,root)
      CALL bcast_real(deltaz,1,root)
!
! ... Sampling namelist ..................................................
!
      IF (mpime == root) READ(punit, sampling)
      CALL bcast_integer(isamp,1,root)
      CALL bcast_integer(number_of_probes,1,root)
      CALL bcast_logical(assign_index,1,root)
      CALL bcast_character(probe_file,80,root)
!
! ... Masspart namelist ................................................
!
      IF (mpime == root) READ(punit, masspart)
      CALL bcast_integer(imassn,1,root)
      CALL bcast_integer(number_of_boxes,1,root)
      CALL bcast_character(boxes_file,80,root)
!
! ... Massflux namelist .................................................
!
      IF (mpime == root) READ(punit, massflux)
      CALL bcast_integer(ifluxn,1,root)
      CALL bcast_integer(number_of_planes,1,root)
      CALL bcast_character(planes_file,80,root)
!
! ... Massgsedim namelist ...............................................
!
      IF (mpime == root) READ(punit, massgsedim)
      CALL bcast_integer(iground,1,root)
      CALL bcast_real(thickness,1,root)
!
! ... Animation namelist ...............................................
!
      IF (mpime == root) READ(punit, animation)
      CALL bcast_integer(variable_n,1,root)
!
! ... Post-processing namelist .........................................
!
      IF (mpime == root) READ(punit, post_processing)
      CALL bcast_integer(field_n,1,root)
!
      RETURN
      END SUBROUTINE postin
!----------------------------------------------------------------------
      END MODULE postp_input
!----------------------------------------------------------------------
