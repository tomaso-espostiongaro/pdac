!----------------------------------------------------------------------
      MODULE process_outp
!----------------------------------------------------------------------
      USE control_flags, ONLY: job_type, formatted_output
      USE dimensions, ONLY: ntot, nz
      USE io_files, ONLY: logunit, tempunit, testunit
      USE kinds
      USE postp_output, ONLY: write_topo2d, write_map
      USE postp_output, ONLY: write_fields, write_mean_field
      USE postp_output, ONLY: write_axis_average
      USE postp_output, ONLY: write_plume_profile
      USE postp_output, ONLY: read_implicit_profile, read_output, read_old_output
      USE postp_output, ONLY: first_out, last_out, incr_out
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
      USE domain_mapping, ONLY: ncdom, ncint
      USE grid, ONLY: z
      USE sample_points, ONLY: isamp, sample
      USE sample_points, ONLY: column_axis_average, allocate_column_probes, column_axis_std
      USE sample_points, ONLY: compute_plume_profile, allocate_plume_variables
      USE mass_partition, ONLY: imassn, massn
      USE mass_orthoflux, ONLY: ifluxn, fluxn
      USE mass_ground, ONLY: iground, massgs
      USE parallel, ONLY: mpime, root
      USE postp_variables, ONLY: allocate_main_fields, allocate_derived_fields, tg
      USE postp_variables, ONLY: compute_derived_fields, allocate_derived_fields_av
      USE postp_variables, ONLY: pd, tg, ts, lepst, rhom, um, vm, wm, mvm, tm
      USE postp_variables, ONLY: rlk, rgp, vf, rhog, p, tg, xgc, mg
      USE gas_constants, ONLY: rgas
      USE postp_variables, ONLY: rhom_gav, um_gav, vm_gav, wm_gav, tm_gav
!
      IMPLICIT NONE
      INTEGER :: tn, cnt, itn
      INTEGER :: ijk, i, j, k, ig, is, n, imesh
      REAL*8  :: md, md1
!
      LOGICAL :: lform, ex
      CHARACTER(LEN = 14) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
!
      lform = formatted_output
!
! ... Loop counters...
!
      cnt = (last_out - first_out) / incr_out + 1
      md = 1.D0 / cnt
      IF (cnt > 1) THEN
        md1 = 1.D0 / (cnt-1)
      ELSE
        md1 = 0.D0
      END IF
!
! ... Each processor allocates the main fields
! ... and the derived fields in the subdomain
!
      CALL allocate_main_fields(ncdom)
      CALL allocate_derived_fields(ncdom)
      IF (imnfld > 0) CALL allocate_derived_fields_av(ncdom) 
      IF (isamp > 1) CALL allocate_column_probes(nz,cnt)
      IF (isamp > 2) CALL allocate_plume_variables(nz)
!
      itn = 0
      DO tn = first_out, last_out, incr_out
        itn = itn + 1
        !
        IF (mpime == root) &
          WRITE(logunit,fmt="(/,'* Starting post-processing ',I5,' * ')" ) tn
        !
        ! ... Read PDAC output file
        !
        CALL read_output ( tn )
        !CALL read_old_output ( tn )
        !
        ! ... Derived fields are computed as a function of
        ! ... primary fields and other derived fields
        !
        IF (mpime == root) &
          WRITE(logunit,fmt="('  Computing derived fields' )")
        CALL compute_derived_fields
        !
        ! ... Write out fields of interest
        !
        IF (iflds > 0) CALL write_fields(tn)
        !
        ! ... Compute time-averaged fields
        !
        IF (imnfld > 0) CALL compute_time_averages
        !
        ! ... Print the map of any interesting variable above ground
        !
        IF (imap > 0) THEN
               CALL write_topo2d
               CALL write_map(tn,lepst,'le')
!               CALL write_map(tn,tg,'tg')
!               CALL write_map(tn,pd,'pd')
        END IF
        IF (isamp > 0)  CALL sample
        IF (isamp > 1)  CALL column_axis_average(itn)
        IF (isamp > 2)  THEN
          CALL compute_plume_profile
          CALL write_plume_profile(tn)
        END IF
        IF (imassn  > 0)  CALL massn
        IF (ifluxn  > 0)  CALL fluxn
        IF (iground > 0)  CALL massgs(tn)
        !
      END DO
!
! ... Compute the standard deviations and
! ... print the time-averaged column probes
!
      IF (isamp > 1) THEN
        CALL column_axis_std(first_out, last_out, incr_out, cnt, md, md1)
        CALL write_axis_average(md,md1)
      END IF
!
! ... Print the time-averaged field 
!
      IF (imnfld > 0) THEN
        rhom_gav = rhom_gav * md
        um_gav = um_gav * md
        IF (job_type == '3D') vm_gav = vm_gav * md
        wm_gav = wm_gav * md
        tm_gav = tm_gav * md
        CALL write_mean_field(last_out)
      END IF
!
      RETURN
      END SUBROUTINE process
!-----------------------------------------------------------------------
      END MODULE process_outp
!----------------------------------------------------------------------
