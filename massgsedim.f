!----------------------------------------------------------------------
      MODULE mass_ground
!----------------------------------------------------------------------
      USE dimensions
      USE filter_outp, ONLY: first_out, last_out, incr_out
      USE filter_outp, ONLY: read_array, write_array
      USE kinds
      USE control_flags, ONLY: job_type
      USE io_files, ONLY: logunit
      USE output_dump, ONLY: formatted_output
!
      IMPLICIT NONE
!
! ... main fields
!
      REAL, ALLOCATABLE, DIMENSION(:)   :: p, ug, vg, wg, tg, eptemp
      REAL, ALLOCATABLE, DIMENSION(:,:) :: xgc
      REAL, ALLOCATABLE, DIMENSION(:,:) :: eps, us, vs, ws, ts
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: m_map, t_map
!
      REAL*8   :: time
      REAL*4   :: stime
      INTEGER  :: iground
      REAL*8   :: thickness
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_main_fields(dime)
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: dime

      ALLOCATE(p(dime), ug(dime), vg(dime), wg(dime), tg(dime))
      ALLOCATE(eptemp(dime))
      ALLOCATE(eps(dime,nsolid), us(dime,nsolid), vs(dime,nsolid), &
                ws(dime,nsolid), ts(dime,nsolid))
      ALLOCATE(xgc(dime,ngas))
      ALLOCATE(m_map(nx,ny))
      ALLOCATE(t_map(nx,ny))
      m_map = 0.D0
      t_map = 0.D0

      RETURN
      END SUBROUTINE allocate_main_fields
!----------------------------------------------------------------------
      SUBROUTINE read_output( tn )
!
! ... Read THe PDAC Output files
!
      USE io_files, ONLY: outpunit
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: tn
      CHARACTER(LEN = 11) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
      CHARACTER(LEN = 2 ) :: lettera2
      LOGICAL :: lform
!
      INTEGER :: ig, is, i, k, ijk
!
      filnam='output.'//lettera(tn)
!
! ... OLD OUTPUT name
!
!      filnam='OUTPUT.'//lettera(tn)

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
!
! ... OLD OUTPUT format had volume fractions here ...
!
!      DO is = 1, nsolid
!
!        IF( lform ) READ(outpunit,'(///)')
!        CALL read_array( outpunit, eps(:,is), lform )  ! solid_bulk_density
!
!      END DO
!
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

      ! ... OLD format
      !DO ig=1,1
!      
! ... NEW output contains all gas species in
!
      DO ig=1,ngas
        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, xgc(:,ig), lform )  ! gc_molar_fraction
      END DO

      DO is = 1, nsolid
!
! ... NEW output format has volume fractions here ...
!
        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, eps(:,is), lform )  ! solid_bulk_density
!
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
! ... OLD OUTPUT were written starting from the domain top
!
!      DO is = 1, nsolid
!      eptemp = 0.D0
!      DO k = 1, nz
!        DO i = 1, nx
!          ijk = i + (k-1) * nx
!          eptemp(ijk) = eps(i + (nz-k)*nx,is)
!          END DO
!        END DO
!        eps(:,is) = eptemp(:)
!      END DO
!
      RETURN
      END SUBROUTINE read_output
!-----------------------------------------------------------------------
      SUBROUTINE massgs
! 
      USE dimensions
      USE grid, ONLY: dx, dy, dz, r, xb, yb, zb, z
      USE immersed_boundaries, ONLY: immb
      USE io_files, ONLY: tempunit
      USE particles_constants, ONLY: rl
      USE filter_outp, ONLY: improfile

      IMPLICIT NONE
!
      INTEGER :: i, j, k, ijk
      INTEGER :: igroud, tn
      REAL*8  :: mmap, quota
      LOGICAL :: lform
      CHARACTER(LEN = 14) :: filnam1, filnam2
      CHARACTER(LEN = 4 ) :: lettera
      REAL*8 :: volume
      REAL*8, ALLOCATABLE :: vf(:)
      REAL*8, ALLOCATABLE :: tmap(:)

      IF (job_type == '2D') RETURN

      ALLOCATE(tmap(nsolid))
      tmap = 0.D0
!
      ALLOCATE(vf(ntot))
      vf(:) = 1.D0
      IF (immb >= 1) THEN
            OPEN(tempunit,FILE='vf.dat',STATUS='OLD',FORM='UNFORMATTED')
            READ(tempunit) (vf(ijk), ijk=1,ntot)
            CLOSE(tempunit)
      END IF
!
      lform = formatted_output

      CALL allocate_main_fields(ntot)
!
      DO tn = first_out, last_out, incr_out

        filnam1='m_map.'//lettera(tn)
        filnam2='t_map.'//lettera(tn)
        WRITE(logunit,fmt="(/,'* Starting post-processing ',I5,' * ')" ) tn

        ! ... Read PDAC output file
        !
        CALL read_output ( tn )

        ! ... Compute the total solid mass to the ground
        !
        DO i = 1, nx
          DO j = 1, ny
            !
            tmap(:) = 0.D0
            mmap = 0.D0
            DO k = 1, nz
              quota = improfile(i,j,k)
              IF (quota >= 0.D0 .AND. quota<=thickness) THEN 
                  ijk  = i + (j-1) * nx + (k-1) * nx * ny
                  volume = dx(i) * dy(j) * dz(k) * vf(ijk)
                  mmap = mmap + SUM(eps(ijk,:)*rl(:)) * dz(k)
                  tmap(:) = tmap(:) + eps(ijk,:) * dz(k)
                  ! ... Map the value reached at any given position at
                  ! ... given time
              END IF
            END DO
            m_map(i,j) = mmap
            t_map(i,j) = SUM(tmap(:))/0.67D0
            !
          END DO
        END DO
        !
      END DO
!

! ... Print out the map and the new 2D DEM file
!
      OPEN(UNIT=tempunit,FILE=filnam1)
      DO j = 1, ny
          WRITE(tempunit,122) (m_map(i,j), i=1, nx)
      END DO
      CLOSE(tempunit)
      OPEN(UNIT=tempunit,FILE=filnam2)
      DO j = 1, ny
          WRITE(tempunit,122) (t_map(i,j), i=1, nx)
      END DO
      CLOSE(tempunit)

 122  FORMAT(10(1x,G14.6E3))

      RETURN
      END SUBROUTINE massgs
!-----------------------------------------------------------------------
      END MODULE mass_ground
!----------------------------------------------------------------------
