!----------------------------------------------------------------------
      MODULE mass_partition
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
!
      TYPE subdomain 
        INTEGER :: i1
        INTEGER :: i2
        INTEGER :: j1
        INTEGER :: j2
        INTEGER :: k1
        INTEGER :: k2
      END TYPE subdomain 
!
      INTEGER :: number_of_boxes
      CHARACTER(LEN=80) :: boxes_file
      REAL*8   :: time
      REAL*4   :: stime
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
      SUBROUTINE massn
! 
      USE dimensions
      USE grid, ONLY: dx, dy, dz, r, xb, yb, zb
      USE immersed_boundaries, ONLY: immb
      USE io_files, ONLY: tempunit
      USE particles_constants, ONLY: rl

      IMPLICIT NONE
!
      INTEGER :: i1, i2, j1, j2, k1, k2
      INTEGER :: is, ijk, i, j, k, nfil, n, tn, counter
      LOGICAL :: lform
      CHARACTER(LEN = 14) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
      REAL*8 :: volume, pi, twopi, totalmass
      REAL*8, ALLOCATABLE :: vf(:)
      REAL*8, ALLOCATABLE :: smass(:,:)
      TYPE( subdomain ), ALLOCATABLE :: box(:)
      REAL*8 :: x1, x2, y1, y2, z1, z2

      ALLOCATE(vf(ntot))
      vf(:) = 1.D0
      IF (immb >= 1) THEN
              OPEN(tempunit,FILE='vf.dat',STATUS='OLD',FORM='UNFORMATTED')
              READ(tempunit) (vf(ijk), ijk=1,ntot)
              CLOSE(tempunit)
      END IF
!
      pi = 4.D0 * ATAN(1.D0)
      twopi = 2.D0 * pi
!
      lform = formatted_output

      ! ... Allocate boxes.
      ! ... Read the file containing the indexes or coordinates of the boxes
      !
      ALLOCATE(box(number_of_boxes))
      ALLOCATE(smass(number_of_boxes,nsolid))
      box(:)%i1 = 1
      box(:)%i2 = 1
      box(:)%j1 = 1
      box(:)%j2 = 1
      box(:)%k1 = 1
      box(:)%k2 = 1

      CALL allocate_main_fields(ntot)
!
      OPEN(tempunit, FILE=boxes_file, STATUS='OLD')
!
      DO n = 1, number_of_boxes
                IF (job_type == '2D') THEN
                        READ(tempunit,*) x1, x2, z1, z2
                        DO i=1,nx
                          IF (xb(i) < x1) box(n)%i1 = i
                          IF (xb(i) < x2) box(n)%i2 = i
                        END DO
                        DO k=1,nz
                          IF (zb(k) < z1) box(n)%k1 = k
                          IF (zb(k) < z2) box(n)%k2 = k
                        END DO
                        box(n)%j1 = 1
                        box(n)%j2 = 1
                ELSE IF (job_type == '3D') THEN
                        READ(tempunit,*) x1, x2, y1, y2, z1, z2
                        DO i=1,nx
                          IF (xb(i) < x1) box(n)%i1 = i
                          IF (xb(i) < x2) box(n)%i2 = i
                        END DO
                        DO j=1,ny
                          IF (yb(j) < y1) box(n)%j1 = j
                          IF (yb(j) < y2) box(n)%j2 = j
                        END DO
                        DO k=1,nz
                          IF (zb(k) < z1) box(n)%k1 = k
                          IF (zb(k) < z2) box(n)%k2 = k
                        END DO
                END IF
      END DO
      CLOSE(tempunit)
      WRITE(*,*) 'Box limits'
      WRITE(*,*) (box(n)%i1, box(n)%i2, box(n)%j1, box(n)%j2, box(n)%k1, box(n)%k2, &
       n=1,number_of_boxes)
!
      OPEN(tempunit, FILE='masspart.dat', STATUS='UNKNOWN', POSITION='APPEND')
      WRITE(tempunit,*) '****************************'
      DO tn = first_out, last_out, incr_out

        WRITE(logunit,fmt="(/,'* Starting post-processing ',I5,' * ')" ) tn

        ! ... Read PDAC output file
        !
        CALL read_output ( tn )

        ! ... Compute the total solid mass in a box
        !
        smass(:,:) = 0.D0
        totalmass = 0.D0
        DO n = 1, number_of_boxes
          counter = 0
          i1 = box(n)%i1
          i2 = box(n)%i2
          j1 = box(n)%j1
          j2 = box(n)%j2
          k1 = box(n)%k1
          k2 = box(n)%k2
          DO k = k1, k2
            DO j = j1, j2
              DO i = i1, i2
                IF (job_type == '2D') THEN
                  ijk = i + (k-1)*nx
                  volume = r(i) * dx(i) * dz(k) * vf(ijk)
                  volume = twopi * volume
                ELSE IF (job_type == '3D') THEN
                  ijk = i + (j-1)*nx + (k-1)*nx*ny
                  volume = dx(i) * dy(j) * dz(k) * vf(ijk)
                END IF
                IF (i/=1 .AND. i/=nx .AND. k/=1 .AND. k/=nz) THEN
                  IF ((j/=1 .AND. j/=ny) .OR. (job_type == '2D')) THEN
                          counter = counter + 1
                          DO is = 1, nsolid
                            smass(n,is)  = smass(n,is) + eps(ijk,is) * rl(is) * volume
                          END DO
                  END IF
                END IF
              END DO
            END DO
          END DO
          !WRITE(*,*) 'Block ', n, ' Counts: ', counter 
        END DO
        totalmass = SUM(smass(1:3,:))
        WRITE(tempunit,100) &
          time, ((smass(n,is),is=1,nsolid),n=1,number_of_boxes),totalmass
      END DO
      CLOSE(tempunit)
 100  FORMAT(20(G30.15E3))
!
      RETURN
      END SUBROUTINE massn
!-----------------------------------------------------------------------
      END MODULE mass_partition
!----------------------------------------------------------------------
