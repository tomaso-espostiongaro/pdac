!----------------------------------------------------------------------
      MODULE process_outp
!----------------------------------------------------------------------
      USE dimensions
      USE postp_output, ONLY: write_topo2d, write_map
      USE postp_output, ONLY: write_fields, write_mean_field
      USE postp_output, ONLY: write_vertical_profiles
      USE postp_output, ONLY: read_implicit_profile, read_output
      USE postp_output, ONLY: first_out, last_out, incr_out
      USE kinds
      USE control_flags, ONLY: job_type, formatted_output
      USE io_files, ONLY: logunit
!
      IMPLICIT NONE
!
      INTEGER :: act
      INTEGER :: iflds, imnfld, imap
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE process
! 
! ... Compute the derived fields
!
      USE compute_mean_fields, ONLY: compute_time_averages
      USE compute_mean_fields, ONLY: compute_vertical_mean_profiles
      USE dimensions, ONLY: ntot
      USE domain_mapping, ONLY: ncdom, ncint
      USE grid, ONLY: z
      USE io_files, ONLY: tempunit, testunit
      USE sample_points, ONLY: isamp, sample
      USE mass_partition, ONLY: imassn, massn
      USE mass_orthoflux, ONLY: ifluxn, fluxn
      USE mass_ground, ONLY: iground, massgs
      USE parallel, ONLY: mpime, root
      USE postp_variables, ONLY: allocate_main_fields, allocate_derived_fields, tg
      USE postp_variables, ONLY: compute_derived_fields, allocate_derived_fields_av
      USE postp_variables, ONLY: pd, tg, ts, lepst, rhom, um, vm, wm, mvm, tm
      USE postp_variables, ONLY: rhom_av, um_av, vm_av, wm_av, tm_av
!
      IMPLICIT NONE
      INTEGER :: tn, cnt
      INTEGER :: ijk, i, j, k, ig, is, n, imesh
      REAL*8  :: md
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
      IF (imnfld > 0) CALL allocate_derived_fields_av(ncdom) 
!
! ... Loop ...
!
      cnt = (last_out - first_out + 1) / incr_out
      md = 1.D0 / cnt

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
        !
        ! ... Compute time-averaged fields
        !
        IF (imnfld > 0) CALL compute_time_averages
        IF (imnfld > 1) CALL compute_vertical_mean_profiles
        !
        ! ... Print the map of any interesting variable above ground
        !
        IF (imap > 0) THEN
               CALL write_topo2d
               CALL write_map(tn,lepst,'le')
               CALL write_map(tn,pd,'pd')
        END IF
        IF (isamp   > 0)  CALL sample
        IF (imassn  > 0)  CALL massn
        IF (ifluxn  > 0)  CALL fluxn
        IF (iground > 0)  CALL massgs(tn)
        IF (imnfld  > 1)  CALL write_vertical_profiles(tn)
        !
      END DO
!
! ... Print the time-averaged field 
!
      IF (imnfld > 0) THEN
        rhom_av = rhom_av * md
        um_av = um_av * md
        IF (job_type == '3D') vm_av = vm_av * md
        wm_av = wm_av * md
        tm_av = tm_av * md
        CALL write_mean_field(last_out)
      END IF
!
      RETURN
      END SUBROUTINE process
!-----------------------------------------------------------------------
      END MODULE process_outp
!----------------------------------------------------------------------
