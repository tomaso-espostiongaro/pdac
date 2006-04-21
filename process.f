!----------------------------------------------------------------------
      MODULE process_outp
!----------------------------------------------------------------------
      USE dimensions
      USE postp_output, ONLY: write_topo2d, write_map, write_fields
      USE postp_output, ONLY: read_implicit_profile, read_output
      USE postp_output, ONLY: first_out, last_out, incr_out
      USE kinds
      USE control_flags, ONLY: job_type, formatted_output
      USE io_files, ONLY: logunit
!
      IMPLICIT NONE
!
      INTEGER :: act
      INTEGER :: iflds, imap
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE process
! 
! ... Compute the derived fields
!
      USE dimensions, ONLY: ntot
      USE domain_mapping, ONLY: ncdom
      USE grid, ONLY: z
      USE io_files, ONLY: tempunit
      USE sample_points, ONLY: isamp, sample
      USE mass_partition, ONLY: imassn, massn
      USE mass_orthoflux, ONLY: ifluxn, fluxn
      USE mass_ground, ONLY: iground, massgs
      USE parallel, ONLY: mpime, root
      USE postp_variables, ONLY: allocate_main_fields, allocate_derived_fields
      USE postp_variables, ONLY: compute_derived_fields, pd, tg, ts, lepstot
!
      IMPLICIT NONE
      INTEGER :: tn, nv, cnt
      INTEGER :: ijk, i, j, k, ig, is, n
!
      LOGICAL :: lform, ex
      CHARACTER(LEN = 14) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
!
      lform = formatted_output
!
! ... Each processor allocates the main fields
! ... and the derived fields in the subdomain
!
      CALL allocate_main_fields(ncdom)
      CALL allocate_derived_fields(ncdom)
!
      DO tn = first_out, last_out, incr_out
        !
        IF (mpime == root) &
          WRITE(logunit,fmt="(/,'* Starting post-processing ',I5,' * ')" ) tn
        !
        ! ... Read PDAC output file
        !
        CALL read_output ( tn )
        !
        ! ... Derived fields are computed as a function of
        ! ... primary fields and other derived fields
        !
        CALL compute_derived_fields
        !
        ! ... Write out fields of interest
        !
        IF (iflds > 0) CALL write_fields(tn)

        ! ... Print the map of any interesting variable above ground
        !
        IF (imap > 0) THEN
               CALL write_topo2d
               CALL write_map(tn,pd,'pd')
               CALL write_map(tn,tg,'tg')
               CALL write_map(tn,ts(:,1),'t1')
               CALL write_map(tn,ts(:,2),'t2')
               CALL write_map(tn,lepstot,'lt')
        END IF
        IF (isamp   > 0)  CALL sample
        IF (imassn  > 0)  CALL massn
        IF (ifluxn  > 0)  CALL fluxn
        IF (iground > 0)  CALL massgs(tn)
        !
      END DO
!
      RETURN
      END SUBROUTINE process
!-----------------------------------------------------------------------
      END MODULE process_outp
!----------------------------------------------------------------------
