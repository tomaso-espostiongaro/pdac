!----------------------------------------------------------------------
      MODULE process_outp
!----------------------------------------------------------------------
!
      USE dimensions
      USE io_serial, ONLY: read_array, write_array
      USE io_serial, ONLY: write_topo2d, write_map
      USE io_serial, ONLY: read_implicit_profile
      USE io_serial, ONLY: first_out, last_out, incr_out
      USE kinds
      USE control_flags, ONLY: job_type
      USE io_files, ONLY: logunit
      USE output_dump, ONLY: formatted_output
      USE postp_variables
!
      IMPLICIT NONE
!
      INTEGER :: act
      INTEGER :: iflds, imap
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE read_output( tn )
!
! ... Read THe PDAC Output files
!
      USE io_files, ONLY: outpunit
!
      IMPLICIT NONE
      SAVE
!
      INTEGER, INTENT(IN) :: tn
      CHARACTER(LEN = 11) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
      CHARACTER(LEN = 2 ) :: lettera2
      LOGICAL :: lform
!
      INTEGER :: ig, is, ii
      REAL*4   :: stime
!
      filnam='output.'//lettera(tn)

      lform = formatted_output
      IF (lform) THEN
        OPEN(UNIT=outpunit, FILE=filnam, STATUS='OLD')
        READ(outpunit,'(1x,///,1x,"@@@ TIME = ",g11.4)') time
      ELSE 
        OPEN(UNIT=outpunit,FORM='UNFORMATTED',FILE=filnam)
        READ(outpunit) stime
        time = stime
      END IF

      IF( lform ) READ(outpunit,'(///)')
      CALL read_array( outpunit, p, lform )  ! gas_pressure

      IF (job_type == '2D') THEN

        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, ug, lform ) ! gas_velocity_r
        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, wg, lform ) ! gas_velocity_z

      ELSE IF (job_type == '3D') THEN

        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, ug, lform ) ! gas_velocity_x
        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, vg, lform ) ! gas_velocity_y
        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, wg, lform ) ! gas_velocity_z

      ELSE
        CALL error('outp_','Unknown job type',1)
      END IF

      IF( lform ) READ(outpunit,'(///)')
      CALL read_array( outpunit, tg, lform )  ! gas_temperature

      DO ig=1,ngas
        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, xgc(:,ig), lform )  ! gc_molar_fraction
      END DO

      DO is = 1, nsolid

        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, eps(:,is), lform )  ! solid_bulk_density

        IF (job_type == '2D') THEN

        IF( lform ) READ(outpunit,'(///)')
          CALL read_array( outpunit, us(:,is), lform )  ! solid_velocity_r
        IF( lform ) READ(outpunit,'(///)')
          CALL read_array( outpunit, ws(:,is), lform )  ! solid_velocity_z

        ELSE IF (job_type == '3D') THEN

        IF( lform ) READ(outpunit,'(///)')
          CALL read_array( outpunit, us(:,is), lform )  ! solid_velocity_x
        IF( lform ) READ(outpunit,'(///)')
          CALL read_array( outpunit, vs(:,is), lform )  ! solid_velocity_y
        IF( lform ) READ(outpunit,'(///)')
          CALL read_array( outpunit, ws(:,is), lform )  ! solid_velocity_z

        END IF

        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, ts(:,is), lform )  ! solid_temperature

      END DO

      CLOSE (outpunit)
!
      RETURN
      END SUBROUTINE read_output
!-----------------------------------------------------------------------
      SUBROUTINE logtg_fields (tn)
      USE io_files, ONLY: tempunit
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: tn
      CHARACTER(LEN = 14) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
      LOGICAL :: lform
!
      filnam = 'log10epst.'//lettera(tn)
      IF (lform) THEN
        OPEN(tempunit,FILE=filnam)
      ELSE
        OPEN(tempunit,FILE=filnam, FORM='UNFORMATTED')
      END IF
      CALL write_array( tempunit, lepstot, lform )
      CLOSE(tempunit)
      !
      filnam = 'tg.'//lettera(tn)
      IF (lform) THEN
        OPEN(tempunit,FILE=filnam)
      ELSE
        OPEN(tempunit,FILE=filnam, FORM='UNFORMATTED')
      END IF
      CALL write_array( tempunit, tg, lform )
      CLOSE(tempunit)
!
      RETURN
      END SUBROUTINE logtg_fields
!-----------------------------------------------------------------------
      SUBROUTINE process
! 
! ... Compute the derived fields
!
      USE dimensions, ONLY: ntot
      USE grid, ONLY: z
      USE io_files, ONLY: tempunit
      USE sample_points, ONLY: isamp, sample
      USE mass_partition, ONLY: imassn, massn
      USE mass_orthoflux, ONLY: ifluxn, fluxn
      USE mass_ground, ONLY: iground, massgs

      IMPLICIT NONE
      INTEGER :: tn, nv, cnt
      INTEGER :: nfil
      INTEGER :: ijk, i, j, k, ig, is, n
!
      LOGICAL :: lform, ex
      CHARACTER(LEN = 14) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
!
      lform = formatted_output
!
      CALL allocate_main_fields(ntot)
      CALL allocate_derived_fields(ntot)
!
      DO tn = first_out, last_out, incr_out
        !
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
        IF (iflds > 0) CALL logtg_fields(tn)

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
