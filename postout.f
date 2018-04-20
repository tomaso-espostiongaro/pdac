!-----------------------------------------------------------------------
      MODULE postp_output
!-----------------------------------------------------------------------
      USE dimensions
      USE kinds
      IMPLICIT NONE
!
      REAL*8, ALLOCATABLE :: improfile_2d(:,:)
      REAL*8, ALLOCATABLE :: improfile_3d(:,:,:)
      REAL*8, ALLOCATABLE, DIMENSION(:) :: map1d
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: map2d
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: topo2d

      INTEGER :: first_out, last_out, incr_out
      ! ... crop
      INTEGER :: iminc, imaxc, jminc, jmaxc, kminc, kmaxc
      LOGICAL :: print_log, print_tg, print_mn, print_cm, print_pd, print_mnn , print_rhom
      REAL*8  :: deltaz1, deltaz2
!
      PUBLIC
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE read_old_output( nf )
      USE control_flags, ONLY: job_type, formatted_output
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nsolid, ngas
      USE gas_constants, ONLY: gas_type
      USE kinds
      USE io_files, ONLY: filnam, outpunit, logunit
      USE io_parallel, ONLY: read_array
      USE postp_variables, ONLY: time, p,ug,vg,wg,tg,xgc,us,vs,ws,ts,eps
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: nf
      CHARACTER( LEN =  4 ) :: lettera
      LOGICAL :: lform
      INTEGER :: is
      REAL*4 :: time4
      REAL*8, ALLOCATABLE :: otmp(:)
!
      filnam='OUTPUT.'//lettera(nf)
      lform = .TRUE.

      IF( mpime == root ) THEN

        OPEN(UNIT=outpunit,FILE=filnam)
        READ(outpunit,155) time
 155    FORMAT(1x,///,12x,g11.4)

        WRITE(logunit,fmt="('  from process: reading file ',A20)") filnam
        WRITE(logunit,*) 'time = ', time
 
      END IF
!
      CALL bcast_real(time, 1, root)
!
      IF( mpime == root ) READ(outpunit,122)
      p = 0.D0
      CALL read_array( outpunit, p, sgl, lform )  ! gas_pressure
      p = p / 10.D0

      ALLOCATE( otmp( SIZE( eps, 1 ) ) )
      DO is = 1, nsolid

        IF( mpime == root ) READ(outpunit,122)
        otmp = 0.D0
        CALL read_array( outpunit, otmp, sgl, lform )  ! solid_volume_fraction
        eps(:,is) = otmp

      END DO
      DEALLOCATE( otmp )

      IF( mpime == root ) READ(outpunit,122)
      wg = 0.D0
      CALL read_array( outpunit, wg, sgl, lform ) ! gas_velocity_z
      wg = wg / 100.D0

      IF( mpime == root ) READ(outpunit,122)
      ug = 0.D0
      CALL read_array( outpunit, ug, sgl, lform ) ! gas_velocity_r
      ug = ug / 100.D0

      IF( mpime == root ) READ(outpunit,122)
      tg = 0.D0
      CALL read_array( outpunit, tg, sgl, lform )  ! gas_temperature
!
      ALLOCATE( otmp( SIZE( xgc, 1 ) ) )
      otmp = 0.D0
      IF( mpime == root ) READ(outpunit,122)
      CALL read_array( outpunit, otmp, sgl, lform )  ! gc_molar_fraction
      xgc(:,1) = otmp
      xgc(:,2) = 1.D0 - otmp
      DEALLOCATE( otmp )
!
      DO is = 1, nsolid

        IF( mpime == root ) READ(outpunit,122)
        ws(:,is) = 0.D0
        CALL read_array( outpunit, ws(:,is), sgl, lform )  ! solid_velocity_z
        ws(:,is)  = ws(:,is) / 100.D0

        IF( mpime == root ) READ(outpunit,122)
        us(:,is) = 0.D0
        CALL read_array( outpunit, us(:,is), sgl, lform )  ! solid_velocity_r
        us(:,is)  = us(:,is) / 100.D0

        IF( mpime == root ) READ(outpunit,122)
        ts(:,is) = 0.D0
        CALL read_array( outpunit, ts(:,is), sgl, lform )  ! solid_temperature

      END DO

      IF( mpime == root ) THEN
        CLOSE (outpunit)
      END IF
!
 122  FORMAT(1x,//,6x,/)

      RETURN
      END SUBROUTINE read_old_output
!----------------------------------------------------------------------
      SUBROUTINE read_output( nf )
      USE control_flags, ONLY: job_type, formatted_output
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nsolid, ngas
      USE gas_constants, ONLY: gas_type
      USE kinds
      USE io_files, ONLY: filnam, outpunit, logunit, testunit
      USE io_parallel, ONLY: read_array
      USE postp_variables, ONLY: time, p,ug,vg,wg,tg,xgc,us,vs,ws,ts,eps
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: nf
      CHARACTER( LEN =  4 ) :: lettera
      LOGICAL :: lform
      INTEGER :: ig,is
      REAL*4 :: time4
      REAL*8, ALLOCATABLE :: otmp(:)
!
      filnam='output.'//lettera(nf)
      lform = formatted_output

      IF( mpime == root ) THEN

        IF (lform) THEN
          OPEN(UNIT=outpunit,FILE=filnam)
          READ(outpunit,155) time
        ELSE 
          OPEN(UNIT=outpunit,FORM='UNFORMATTED',FILE=filnam)
          READ(outpunit) time4
          WRITE(logunit,*) 'Reading time: ', time4
          time = REAL(time4,dbl)
        END IF
 155    FORMAT(1x,///,12x,g11.4)

        WRITE(logunit,fmt="('  from process: reading file ',A20)") filnam
        WRITE(logunit,*) 'time = ', time
 
      END IF
!
      CALL bcast_real(time, 1, root)
!
      IF( lform .AND. mpime == root ) READ(outpunit,122)
      p = 0.D0
      CALL read_array( outpunit, p, sgl, lform )  ! gas_pressure

      IF (job_type == JOB_TYPE_2D) THEN

        IF( lform .AND. mpime == root ) READ(outpunit,122)
        ug = 0.D0
        CALL read_array( outpunit, ug, sgl, lform ) ! gas_velocity_r

        IF( lform .AND. mpime == root ) READ(outpunit,122)
        wg = 0.D0
        CALL read_array( outpunit, wg, sgl, lform ) ! gas_velocity_z

      ELSE IF (job_type == JOB_TYPE_3D) THEN

        IF( lform .AND. mpime == root ) READ(outpunit,122)
        ug = 0.D0
        CALL read_array( outpunit, ug, sgl, lform ) ! gas_velocity_x

        IF( lform .AND. mpime == root ) READ(outpunit,122)
        vg = 0.D0
        CALL read_array( outpunit, vg, sgl, lform ) ! gas_velocity_y

        IF( lform .AND. mpime == root ) READ(outpunit,122)
        wg = 0.D0
        CALL read_array( outpunit, wg, sgl, lform ) ! gas_velocity_z

      ELSE
        CALL error('outp_','Unknown job type',1)
      END IF

      IF( lform .AND. mpime == root ) READ(outpunit,122)
      tg = 0.D0
      CALL read_array( outpunit, tg, sgl, lform )  ! gas_temperature
!
      ALLOCATE( otmp( SIZE( xgc, 1 ) ) )
      DO ig=1,ngas
          IF( lform .AND. mpime == root ) READ(outpunit,122)
          otmp = 0.D0
          CALL read_array( outpunit, otmp, sgl, lform )  ! gc_molar_fraction
          xgc(:,ig) = otmp
      END DO
      DEALLOCATE( otmp )
!
      ALLOCATE( otmp( SIZE( eps, 1 ) ) )

      DO is = 1, nsolid

        IF( lform .AND. mpime == root ) READ(outpunit,122)
        otmp = 0.D0
        CALL read_array( outpunit, otmp, sgl, lform )  ! solid_volume_fraction
        eps(:,is) = otmp

        IF (job_type == JOB_TYPE_2D) THEN

          IF( lform .AND. mpime == root ) READ(outpunit,122)
          us(:,is) = 0.D0
          CALL read_array( outpunit, us(:,is), sgl, lform )  ! solid_velocity_r

          IF( lform .AND. mpime == root ) READ(outpunit,122)
          ws(:,is) = 0.D0
          CALL read_array( outpunit, ws(:,is), sgl, lform )  ! solid_velocity_z

        ELSE IF (job_type == JOB_TYPE_3D) THEN

          IF( lform .AND. mpime == root ) READ(outpunit,122)
          us(:,is) = 0.D0
          CALL read_array( outpunit, us(:,is), sgl, lform )  ! solid_velocity_x

          IF( lform .AND. mpime == root ) READ(outpunit,122)
          vs(:,is) = 0.D0
          CALL read_array( outpunit, vs(:,is), sgl, lform )  ! solid_velocity_y

          IF( lform .AND. mpime == root ) READ(outpunit,122)
          ws(:,is) = 0.D0
          CALL read_array( outpunit, ws(:,is), sgl, lform )  ! solid_velocity_z

        END IF

        IF( lform .AND. mpime == root ) READ(outpunit,122)
        ts(:,is) = 0.D0
        CALL read_array( outpunit, ts(:,is), sgl, lform )  ! solid_temperature

      END DO

      DEALLOCATE( otmp )

      IF( mpime == root ) THEN
        CLOSE (outpunit)
      END IF
!
 122  FORMAT(1x,//,6x,/)

      RETURN
!
      END SUBROUTINE read_output
!----------------------------------------------------------------------
      SUBROUTINE read_implicit_profile
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE grid, ONLY: x, y, z, xb, yb, zb
      USE io_files, ONLY: tempunit
      USE parallel, ONLY: mpime, root

      !  This subroutine reads the implicit profile
      !  and the mesh file computed by PDAC when a 3D
      !  volcano topography is imported from a DEM file

      INTEGER :: i,j,k
      IF (job_type == JOB_TYPE_2D) THEN
        ALLOCATE(improfile_2d(nx,nz))
      ELSE IF (job_type == JOB_TYPE_3D) THEN
        ALLOCATE(improfile_3d(nx,ny,nz))
      END IF
!
! ... Read the georeferenced mesh
!
      IF (mpime == root) THEN
        OPEN(tempunit,FILE='mesh.dat',STATUS='OLD')
        READ(tempunit,*)
        READ(tempunit,*)
        READ(tempunit,*) (x(i), i=1,nx)
        READ(tempunit,*)
        READ(tempunit,*) (xb(i), i=1,nx)
        READ(tempunit,*)
        READ(tempunit,*) (y(j), j=1,ny)
        READ(tempunit,*)
        READ(tempunit,*) (yb(j), j=1,ny)
        READ(tempunit,*)
        READ(tempunit,*) (z(k), k=1,nz)
        READ(tempunit,*)
        READ(tempunit,*) (zb(k), k=1,nz)
        CLOSE(tempunit)
      END IF
!
      CALL bcast_real(x,nx,root)
      CALL bcast_real(xb,nx,root)
      CALL bcast_real(y,ny,root)
      CALL bcast_real(yb,ny,root)
      CALL bcast_real(z,nz,root)
      CALL bcast_real(zb,nz,root)
!      
! ... Read the Implicit Profile
!
      IF (mpime == root) THEN
        OPEN(tempunit,                    &
             FILE='improfile.dat',        &
             STATUS='OLD',                &
             FORM='UNFORMATTED')
        IF (job_type == JOB_TYPE_2D) THEN
          !DO k=1,nz
              !DO i=1,nx
                !READ(tempunit,*) improfile_2d(i,k)
              !END DO
          !END DO
          READ(tempunit) ((improfile_2d(i,k), k=1,nz), i=1,nx)
        ELSE IF (job_type == JOB_TYPE_3D) THEN
          !DO k=1,nz
            !DO j=1,ny
              !DO i=1,nx
                !READ(tempunit,*) improfile_3d(i,j,k)
              !END DO
            !END DO
          !END DO
          READ(tempunit) (((improfile_3d(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        END IF
        CLOSE(tempunit)
      END IF
!
      IF (job_type == JOB_TYPE_2D) THEN
        CALL bcast_real(improfile_2d,ntot,root)
      ELSE IF (job_type == JOB_TYPE_3D) THEN
        CALL bcast_real(improfile_3d,ntot,root)
      END IF
!
 100  FORMAT(5(F20.6))
      RETURN
!
      END SUBROUTINE read_implicit_profile
!-----------------------------------------------------------------------
      SUBROUTINE write_topo2d
! ... Read the implicit profile and write out a DEM file
! ... Notice that if the immersed boundaries are not used,
! ... the DEM is defined by the vertical mesh resolution
!
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nx, ny, nz
      USE grid, ONLY: z, zb, dz
      USE immersed_boundaries, ONLY: immb
      USE io_files, ONLY: tempunit
      USE parallel, ONLY: mpime, root

      IMPLICIT NONE

      REAL*8 :: quota
      INTEGER :: i, j, k, ijk

      IF (job_type == JOB_TYPE_2D) RETURN

      IF (.NOT.ALLOCATED(topo2d)) ALLOCATE(topo2d(nx,ny))
      topo2d = -9999 

      DO i = 2, nx-1, 1
        DO j = 2, ny-1, 1
          !
          search: DO k = 2, nz
            quota = improfile_3d(i,j,k)
            IF (immb > 0 .AND. quota >= 0.D0 .AND. topo2d(i,j) == -9999) THEN
                topo2d(i,j) = z(k) - quota
                EXIT search
            ELSE IF (immb == 0 .AND. quota > -0.5D0*dz(k) .AND. topo2d(i,j) == -9999) THEN
                topo2d(i,j) = zb(k-1)
                EXIT search
            END IF
          END DO search
          !
        END DO
      END DO
      topo2d(1,:) = topo2d(2,:)
      topo2d(nx,:) = topo2d(nx-1,:)
      topo2d(:,1) = topo2d(:,2)
      topo2d(:,ny) = topo2d(:,ny-1)

!
! ... Print out the new 2D DEM file
!
      IF (mpime == root) THEN
        OPEN(UNIT=tempunit,FILE='topo2d.ascii')
        DO j = 1, ny
            WRITE(tempunit,122) (topo2d(i,j), i=1, nx)
        END DO
        CLOSE(tempunit)
      END IF

 122  FORMAT(10(1x,G14.6E3))

      RETURN
      END SUBROUTINE write_topo2d
!-----------------------------------------------------------------------
      SUBROUTINE write_map(nf,array,labl)
      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:) :: array
      INTEGER, INTENT(IN) :: nf
      CHARACTER(LEN=2), INTENT(IN) :: labl
      
      CALL write_map_delta(nf,array,labl,deltaz1)
      CALL write_map_delta(nf,array,labl,deltaz2)

      RETURN
      END SUBROUTINE write_map
!-----------------------------------------------------------------------
      SUBROUTINE write_map_delta(nf,array,labl,deltaz)

      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nx, ny, nz
      USE domain_decomposition, ONLY: cell_owner, cell_g2l
      USE grid, ONLY: z, x, y, flag
      USE io_files, ONLY: tempunit, testunit
      USE parallel, ONLY: mpime, root
      USE set_indexes, ONLY: first_subscr, ijkm
      IMPLICIT NONE

      REAL*8, INTENT(IN) :: deltaz
      REAL*8, INTENT(IN), DIMENSION(:) :: array
      INTEGER, INTENT(IN) :: nf
      CHARACTER(LEN=2), INTENT(IN) :: labl
      REAL*8 :: alpha, map, quota, quotam, value
      INTEGER :: i, j, k, ijk, imesh
      CHARACTER( LEN = 4 ) :: lettera
      CHARACTER( LEN = 20 ) :: filnam
      LOGICAL :: nointerp

      IF (job_type == JOB_TYPE_2D) THEN
!
        ALLOCATE(map1d(nx))
        map1d = 0.D0
      
        filnam='map_'//labl//'_'//lettera(INT(deltaz))//'m'//'.'//lettera(nf)
!
! ... NON HA SENSO CHE QUESTA PROCEDURA DI RICERCA VADA RIPETUTA PER
! ... OGNI MAPPA!
        alpha = 0.0D0
        DO i = 2, nx-1
            !
            search1: DO k = 2, nz
              quota = improfile_2d(i,k)
              IF (quota >= deltaz) THEN
                  quotam = improfile_2d(i,k-1)
                  imesh  = i + (k-1) * nx
                  IF (cell_owner(imesh) == mpime) THEN
                    ijk  = cell_g2l(imesh, mpime)
                    !
                    ! ... identify the first neighbours
                    !
                    CALL first_subscr(ijk)
                    alpha = deltaz - quotam
                    alpha = alpha / (z(k) - z(k-1))
                    !
                    IF (quotam > 0.D0) THEN
                      map = alpha* array(ijk) + (1.D0-alpha) * array(ijkm)
                    ELSE
                      map = alpha* array(ijk) + (1.D0-alpha) * array(ijk)
                    END IF
                    ! ... Map the value reached at any given position at
                    ! ... given time
                    map1d(i) = map
                  END IF
                  EXIT search1
              END IF
            END DO search1
          !
        END DO
!
! ... Assemble the map array
!
        CALL parallel_sum_real(map1d, nx)
        map1d(1) = map1d(2)
        map1d(nx) = map1d(nx-1)
!
! ... Print out the map and the new 2D DEM file
!
        IF (mpime == root) THEN
          OPEN(UNIT=tempunit,FILE=filnam)
          DO i = 2, nx-1
              WRITE (tempunit,121) x(i), map1d(i)
          END DO
          CLOSE(tempunit)
        END IF
        DEALLOCATE(map1d)
!
      ELSE IF (job_type == JOB_TYPE_3D) THEN

        ALLOCATE(map2d(nx,ny))
        map2d = 0.D0
      
        filnam='map_'//labl//'_'//lettera(INT(deltaz))//'m'//'.'//lettera(nf)
!
        alpha = 0.0D0
        DO i = 2, nx-1
        DO j = 2, ny-1
            !
            search2: DO k = 2, nz
              ! ... compute the quota of the cell
              ! ... with respect to the topography
              quota = z(k) - topo2d(i,j)
              imesh  = i + (j-1) * nx + (k-1) * nx * ny
              IF (quota >= deltaz) THEN
                  ! ... compute the quota of the last cell below the threshold
                  ! ... with respect to the topography
                  quotam = z(k-1) - topo2d(i,j)
                  !
                  IF (cell_owner(imesh) == mpime) THEN
                    ijk  = cell_g2l(imesh, mpime)
                    !
                    ! ... identify the first neighbours
                    !
                    CALL first_subscr(ijk)
                    alpha = deltaz - quotam
                    alpha = alpha / (z(k) - z(k-1))
                    !
                    value = array(ijkm)
                    ! NON-STANDARD FORTRAN CALL 'ISNAN' (check Not-A-Number values)
                    !IF (ISNAN(value)) nointerp = .TRUE.
                    IF (quotam > 0.D0) THEN
                      map = alpha* array(ijk) + (1.D0-alpha) * array(ijkm)
                    ELSE
                      map = array(ijk)
                    END IF
                    !IF (nointerp) map = array(ijk)
                    ! ... Map the value reached at any given position at
                    ! ... given time
                    map2d(i,j) = map
                    !IF (ISNAN(map)) CALL error('postout','problem writing map',imesh)
                  END IF
                  EXIT search2
              END IF
            END DO search2
            !
          END DO
        END DO
!
! ... Assemble the map array
!
        CALL parallel_sum_real(map2d, nx*ny)
        map2d(1,:) = map2d(2,:)
        map2d(nx,:) = map2d(nx-1,:)
        map2d(:,1) = map2d(:,2)
        map2d(:,ny) = map2d(:,ny-1)
!
! ... Print out the map and the new 2D DEM file
!
        IF (mpime == root) THEN
          OPEN(UNIT=tempunit,FILE=filnam)
          DO j = 1, ny
              WRITE(tempunit,122) (map2d(i,j), i=1, nx)
          END DO
          CLOSE(tempunit)
        END IF
        DEALLOCATE(map2d)
!
      END IF
!
 121  FORMAT(G12.4,G14.6E3)
 122  FORMAT(10(1x,G14.6E3))
!
      RETURN
      END SUBROUTINE write_map_delta
!----------------------------------------------------------------------
      SUBROUTINE write_fields(tn)
      USE control_flags, ONLY: formatted_output
      USE kinds
      USE io_files, ONLY: tempunit
      USE io_parallel, ONLY: write_crop_array
      USE parallel, ONLY: mpime, root
      USE postp_variables, ONLY: epst, lepst, tg, mn, cm, pd, mnn, tm, rhom
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: tn
      CHARACTER(LEN = 14) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
      LOGICAL :: lform
!
      lform = formatted_output
!
      IF (print_log) THEN
        IF (mpime == root) THEN
          filnam = 'lepst.'//lettera(tn)
          IF (lform) THEN
            OPEN(tempunit,FILE=filnam)
          ELSE
            OPEN(tempunit,FILE=filnam, FORM='UNFORMATTED')
          END IF
        END IF
        !
        CALL write_crop_array( tempunit, lepst, sgl, lform, iminc, imaxc, jminc, jmaxc, kminc, kmaxc )
        !
        IF (mpime == root) CLOSE(tempunit)
      END IF
!
      IF (print_tg) THEN
        IF (mpime == root) THEN
          filnam = 'tg.'//lettera(tn)
          IF (lform) THEN
            OPEN(tempunit,FILE=filnam)
          ELSE
            OPEN(tempunit,FILE=filnam, FORM='UNFORMATTED')
          END IF
        END IF
        !
        CALL write_crop_array( tempunit, tm, sgl, lform, iminc, imaxc, jminc, jmaxc, kminc, kmaxc)
        !
        IF (mpime == root) CLOSE(tempunit)
      END IF
!
      IF (print_rhom) THEN
        IF (mpime == root) THEN
          filnam = 'rhom.'//lettera(tn)
          IF (lform) THEN
            OPEN(tempunit,FILE=filnam)
          ELSE
            OPEN(tempunit,FILE=filnam, FORM='UNFORMATTED')
          END IF
        END IF
        !
        CALL write_crop_array( tempunit, rhom, sgl, lform, iminc, imaxc, jminc, jmaxc, kminc, kmaxc)
        !
        IF (mpime == root) CLOSE(tempunit)
      END IF
!
!
      IF (print_mn) THEN
        IF (mpime == root) THEN
          filnam = 'mach.'//lettera(tn)
          IF (lform) THEN
            OPEN(tempunit,FILE=filnam)
          ELSE
            OPEN(tempunit,FILE=filnam, FORM='UNFORMATTED')
          END IF
        END IF
        !
        CALL write_crop_array( tempunit, mn, sgl, lform, iminc, imaxc, jminc, jmaxc, kminc, kmaxc)
        !
        IF (mpime == root) CLOSE(tempunit)
      END IF
!
      IF (print_cm) THEN
        IF (mpime == root) THEN
          filnam = 'cm.'//lettera(tn)
          IF (lform) THEN
            OPEN(tempunit,FILE=filnam)
          ELSE
            OPEN(tempunit,FILE=filnam, FORM='UNFORMATTED')
          END IF
        END IF
        !
        CALL write_crop_array( tempunit, cm, sgl, lform, iminc, imaxc, jminc, jmaxc, kminc, kmaxc )
        !
        IF (mpime == root) CLOSE(tempunit)
      END IF
!
      IF (print_pd) THEN
        IF (mpime == root) THEN
          filnam = 'pd.'//lettera(tn)
          IF (lform) THEN
            OPEN(tempunit,FILE=filnam)
          ELSE
            OPEN(tempunit,FILE=filnam, FORM='UNFORMATTED')
          END IF
        END IF
        !
        CALL write_crop_array( tempunit, pd, sgl, lform, iminc, imaxc, jminc, jmaxc, kminc, kmaxc )
        !
        IF (mpime == root) CLOSE(tempunit)
      END IF
!
      IF (print_mnn) THEN
        IF (mpime == root) THEN
          filnam = 'macn.'//lettera(tn)
          IF (lform) THEN
            OPEN(tempunit,FILE=filnam)
          ELSE
            OPEN(tempunit,FILE=filnam, FORM='UNFORMATTED')
          END IF
        END IF
        !
        CALL write_crop_array( tempunit, mnn, sgl, lform, iminc, imaxc, jminc, jmaxc, kminc, kmaxc )
        !
        IF (mpime == root) CLOSE(tempunit)
      END IF
!
      RETURN
      END SUBROUTINE write_fields
!-----------------------------------------------------------------------
      SUBROUTINE write_mean_field(last_out)
!
! ... Write the mean field of all process fields
!
      USE control_flags, ONLY: formatted_output,  job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE kinds
      USE dimensions, ONLY: nsolid
      USE io_files, ONLY: tempunit
      USE io_parallel, ONLY: write_crop_array
      USE parallel, ONLY: mpime, root
      USE postp_variables, ONLY: rhom_gav, tm_gav, um_gav, vm_gav, wm_gav, time
!
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: last_out

      CHARACTER(LEN = 14) :: filnam
      CHARACTER*4 :: lettera
      LOGICAL :: lform
!
      INTEGER :: ig,is
      REAL*8, ALLOCATABLE :: otmp(:)
!
      filnam='m_outp.'//lettera(last_out)
      lform = formatted_output

      IF (mpime == root ) THEN
        IF (lform) THEN
          OPEN(UNIT=tempunit,FILE=filnam)
          WRITE(tempunit,*) time
        ELSE
          OPEN(UNIT=tempunit,FORM='UNFORMATTED',FILE=filnam)
          WRITE(tempunit) REAL(time,4)
        END IF
      END IF
!
      CALL write_crop_array( tempunit, rhom_gav, sgl, lform, iminc, imaxc, jminc, jmaxc, kminc, kmaxc )
!
      IF (job_type == JOB_TYPE_2D) THEN
        CALL write_crop_array( tempunit, um_gav, sgl, lform, iminc, imaxc, jminc, jmaxc, kminc, kmaxc)
        CALL write_crop_array( tempunit, wm_gav, sgl, lform, iminc, imaxc, jminc, jmaxc, kminc, kmaxc)
      ELSE IF (job_type == JOB_TYPE_3D) THEN
        CALL write_crop_array( tempunit, um_gav, sgl, lform, iminc, imaxc, jminc, jmaxc, kminc, kmaxc)
        CALL write_crop_array( tempunit, vm_gav, sgl, lform, iminc, imaxc, jminc, jmaxc, kminc, kmaxc)
        CALL write_crop_array( tempunit, wm_gav, sgl, lform, iminc, imaxc, jminc, jmaxc, kminc, kmaxc)
      ELSE
        CALL error('outp_','Unknown job type',1)
      END IF
!
      CALL write_crop_array( tempunit, tm_gav, sgl, lform, iminc, imaxc, jminc, jmaxc, kminc, kmaxc)
!
      IF (mpime == root) CLOSE(tempunit)
!
      RETURN

      END SUBROUTINE write_mean_field
!----------------------------------------------------------------------
      SUBROUTINE write_axis_average(md, md1)
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: ngas, nsolid, nz
      USE grid, ONLY: z
      USE io_files, ONLY: tempunit
      USE parallel, ONLY: mpime, root
      USE postp_variables, ONLY: p_av, ug_av, vg_av, wg_av, tg_av, xgc_av, &
         eps_av, us_av, vs_av, ws_av, ts_av, rhom_av, wm_av, pd_av
      USE postp_variables, ONLY: p_sd, ug_sd, vg_sd, wg_sd, tg_sd, xgc_sd, &
         eps_sd, us_sd, vs_sd, ws_sd, ts_sd, rhom_sd, wm_sd, pd_sd
      USE postp_variables, ONLY: p_axis, ug_axis, vg_axis, wg_axis, tg_axis, xgc_axis, &
         eps_axis, us_axis, vs_axis, ws_axis, ts_axis, rhom_axis, wm_axis, pd_axis
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: md, md1
      INTEGER :: nv, k
!
! ... Writing the average arrays on file
! 
      IF (mpime==root) THEN
        OPEN (UNIT=tempunit, FILE='column.dat', STATUS='UNKNOWN')
        DO k = 1, nz
          IF (job_type == JOB_TYPE_2D) THEN
              WRITE(tempunit,100) k, z(k), &
                 p_av(k), p_sd(k), &
                 ug_av(k), ug_sd(k), &
                 wg_av(k), wg_sd(k), &
                 tg_av(k), tg_sd(k), &
                 (xgc_av(k,nv), xgc_sd(k,nv), nv=1, ngas), &
                 (eps_av(k,nv), eps_sd(k,nv), &
                  us_av(k,nv), us_sd(k,nv), &
                  ws_av(k,nv), ws_sd(k,nv), &
                  ts_av(k,nv), ts_sd(k,nv), nv=1, nsolid), &
                 rhom_av(k), rhom_sd(k), &
                 wm_av(k), wm_sd(k), &
                 pd_av(k), pd_sd(k)
           ELSE IF (job_type == JOB_TYPE_3D) THEN
              WRITE(tempunit,100) k, z(k), &
                 p_av(k), p_sd(k), &
                 ug_av(k), ug_sd(k), &
                 vg_av(k), vg_sd(k), &
                 wg_av(k), wg_sd(k), &
                 tg_av(k), tg_sd(k), &
                 (xgc_av(k,nv), xgc_sd(k,nv), nv=1, ngas), &
                 (eps_av(k,nv), eps_sd(k,nv), &
                  us_av(k,nv), us_sd(k,nv), &
                  vs_av(k,nv), vs_sd(k,nv), &
                  ws_av(k,nv), ws_sd(k,nv), &
                  ts_av(k,nv), ts_sd(k,nv), nv=1, nsolid), &
                 rhom_av(k), rhom_sd(k), &
                 wm_av(k), wm_sd(k), &
                 pd_av(k), pd_sd(k)
          END IF
        END DO  
        CLOSE(tempunit)
      END IF
!
 100  FORMAT( I5, F8.2, 100(G14.6E3,1X) )
!
      DEALLOCATE(p_av)
      DEALLOCATE(ug_av)
      DEALLOCATE(vg_av)
      DEALLOCATE(wg_av)
      DEALLOCATE(tg_av)
      DEALLOCATE(xgc_av)
      DEALLOCATE(eps_av)
      DEALLOCATE(us_av)
      DEALLOCATE(vs_av)
      DEALLOCATE(ws_av)
      DEALLOCATE(ts_av)
      DEALLOCATE(rhom_av)
      DEALLOCATE(wm_av)
      DEALLOCATE(pd_av)
!
      DEALLOCATE(p_axis)
      DEALLOCATE(ug_axis)
      DEALLOCATE(vg_axis)
      DEALLOCATE(wg_axis)
      DEALLOCATE(tg_axis)
      DEALLOCATE(xgc_axis)
      DEALLOCATE(eps_axis)
      DEALLOCATE(us_axis)
      DEALLOCATE(vs_axis)
      DEALLOCATE(ws_axis)
      DEALLOCATE(ts_axis)
      DEALLOCATE(rhom_axis)
      DEALLOCATE(wm_axis)
      DEALLOCATE(pd_axis)
!
      DEALLOCATE(p_sd)
      DEALLOCATE(ug_sd)
      DEALLOCATE(vg_sd)
      DEALLOCATE(wg_sd)
      DEALLOCATE(tg_sd)
      DEALLOCATE(xgc_sd)
      DEALLOCATE(eps_sd)
      DEALLOCATE(us_sd)
      DEALLOCATE(vs_sd)
      DEALLOCATE(ws_sd)
      DEALLOCATE(ts_sd)
      DEALLOCATE(rhom_sd)
      DEALLOCATE(wm_sd)
!
      RETURN
!
 299  CALL error ('sample_column','error in opening tempunit', tempunit)
!
      END SUBROUTINE write_axis_average
!----------------------------------------------------------------------
      SUBROUTINE write_plume_profile(nf)

      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nz
      USE grid, ONLY: z
      USE io_files, ONLY: tempunit
      USE parallel, ONLY: mpime, root
      USE postp_variables, ONLY: rhom_z, tm_z, um_z, vm_z, wm_z, time, surface
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nf
      INTEGER :: k
      CHARACTER( LEN = 4 ) :: lettera
      CHARACTER( LEN = 13 ) :: filnam
      REAL*8 :: pi, invpi
!
      pi = 4.D0 * ATAN(1.D0)
      invpi = 1.D0 / pi

!      IF (job_type == JOB_TYPE_2D) RETURN
      
      filnam='zprofile.'//lettera(nf)
!
! ... Print out the profile
!
      IF (mpime == root) THEN
        OPEN(UNIT=tempunit,FILE=filnam)
        DO k = 1, nz
          IF (job_type == JOB_TYPE_3D) THEN
            WRITE(tempunit,122) z(k),rhom_z(k),um_z(k),vm_z(k),wm_z(k),tm_z(k), DSQRT(surface(k)*invpi)
          ELSE 
            WRITE(tempunit,122) z(k),rhom_z(k),um_z(k),wm_z(k),tm_z(k), DSQRT(surface(k)*invpi)
          END IF
        END DO
        CLOSE(tempunit)
      END IF
!
 122  FORMAT(7(1x,G14.6E3))
!
      RETURN
      END SUBROUTINE write_plume_profile
!-----------------------------------------------------------------------
      SUBROUTINE write_AVS_files
      USE control_flags, ONLY: job_type, formatted_output
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions
      USE iotk_module
      USE grid
      USE io_files, ONLY: iuni_scalar, iuni_u, iuni_v, iuni_w
      IMPLICIT NONE
 
      CHARACTER(LEN=15) :: filetype
      INTEGER :: skip_m(3)
      INTEGER :: skip, skip_time, skip_phase, skip_block
      INTEGER :: lbl, ig, is
      INTEGER :: nfields, ndim, veclen, nlx, nly, nlz

      OPEN( UNIT=iuni_scalar, FILE='scalar.fld', STATUS='UNKNOWN')
      OPEN( UNIT=iuni_u, FILE='u.fld', STATUS='UNKNOWN')
      OPEN( UNIT=iuni_v, FILE='v.fld', STATUS='UNKNOWN')
      OPEN( UNIT=iuni_w, FILE='w.fld', STATUS='UNKNOWN')

!
! ... control parameters
!
      IF (job_type == JOB_TYPE_3D) THEN
        ndim = 3
      ELSE IF(job_type == JOB_TYPE_2D) THEN
        ndim = 2
      ELSE
        WRITE(*,*) 'Unknown job_type'
        STOP
      END IF
!
      IF (formatted_output) THEN
        filetype = 'ascii'
      ELSE
        filetype = 'unformatted'
      END IF

! ... Total number of scalar fields
!
      nfields = nphase*2 + ngas
!
      IF (formatted_output) THEN 
        ! ... number of lines for each block
        lbl = nz * (INT(nx*ny / 10) + MIN(1,MOD(nx*ny,10))) + 4
      ELSE
        ! ... number of bytes for each block
        lbl = ntot * 4 + 8
      END IF 

! ... Number of lines in mesh file
!
        nlx = ( nx/5 + MIN(1,MOD(nx,5)) + 1 )
        nly = ( ny/5 + MIN(1,MOD(ny,5)) + 1 )
        nlz = ( nz/5 + MIN(1,MOD(nz,5)) + 1 )

      IF (formatted_output) THEN
        skip_time = 4 + 4
        skip_phase = (lbl * (ndim+2))
        skip_block = lbl
      ELSE
        skip_time = 12
        skip_phase = (lbl * (ndim+2))
        skip_block = lbl
      END IF

!***** S C A L A R S **************************************************
! ... Common Header
!
      WRITE( iuni_scalar, fmt = 1 )
      WRITE( iuni_scalar, fmt = 2 ) ndim
      IF (ndim == 3) THEN
        WRITE( iuni_scalar, fmt = 3 ) nx
        WRITE( iuni_scalar, fmt = 4 ) ny
        WRITE( iuni_scalar, fmt = 5 ) nz
      ELSE IF (ndim == 2) THEN
        WRITE( iuni_scalar, fmt = 3 ) nx
        WRITE( iuni_scalar, fmt = 4 ) nz
      END IF
      WRITE( iuni_scalar, fmt = 6 ) ndim
      WRITE( iuni_scalar, fmt = 7 )
      WRITE( iuni_scalar, fmt = 8 )
!
! ... specify variables and meshes
!
      veclen = nfields
      WRITE( iuni_scalar, fmt = 9 ) veclen
      skip_m(1) = 2
      WRITE( iuni_scalar, fmt = 10 ) skip_m(1)
      IF (ndim == 3) THEN
        skip_m(2) = 2 + 2*nlx
        WRITE( iuni_scalar, fmt = 11 ) skip_m(2)
        skip_m(3) = 2 + 2*nlx + 2*nly
        WRITE( iuni_scalar, fmt = 12 ) skip_m(3)
      ELSE IF (ndim == 2) THEN
        skip_m(3) = 2 + 2*nlx + 2*nly
        WRITE( iuni_scalar, fmt = 11 ) skip_m(3)
      END IF
!
! ... Scalar variables of GAS phase
!
      skip = skip_time
      WRITE( iuni_scalar, fmt = 13 ) 1, filetype, skip
      skip = skip + skip_block * (ndim+1)
      WRITE( iuni_scalar, fmt = 13 ) 2, filetype, skip

! ... Scalar variables of GAS COMPONENTS
!
      DO ig = 1, ngas
        skip = skip_time + skip_phase + (ig-1) * skip_block
        WRITE( iuni_scalar, fmt = 13 ) 2+ig, filetype, skip
      END DO

! ... Scalar variables of SOLID phases
!
      DO is = 1, nsolid
        skip = skip_time + skip_phase * is + ngas * skip_block
        WRITE( iuni_scalar, fmt = 13 ) 2+ngas+(is-1)*2+1, filetype, skip
        skip = skip + skip_block * (ndim+1)
        WRITE( iuni_scalar, fmt = 13 ) 2+ngas+(is)*2, filetype, skip
      END DO

      WRITE( iuni_scalar, fmt = 14 )
!
!***** V E L O C I T Y   X ********************************************
! ... Common Header
!
      WRITE( iuni_u, fmt = 1 )
      WRITE( iuni_u, fmt = 2 ) ndim
      IF (ndim == 3) THEN
        WRITE( iuni_u, fmt = 3 ) nx
        WRITE( iuni_u, fmt = 4 ) ny
        WRITE( iuni_u, fmt = 5 ) nz
      ELSE IF (ndim == 2) THEN
        WRITE( iuni_u, fmt = 3 ) nx
        WRITE( iuni_u, fmt = 4 ) nz
      END IF
      WRITE( iuni_u, fmt = 6 ) ndim
      WRITE( iuni_u, fmt = 7 )
      WRITE( iuni_u, fmt = 8 )
!
! ... specify variables and meshes
!
      veclen = nphase
      WRITE( iuni_u, fmt = 9 ) veclen
      skip_m(1) = 2 + nlx
      WRITE( iuni_u, fmt = 10 ) skip_m(1)
      IF (ndim == 3) THEN
        skip_m(2) = 2 + 2*nlx
        WRITE( iuni_u, fmt = 11 ) skip_m(2)
        skip_m(3) = 2 + 2*nlx + 2*nly
        WRITE( iuni_u, fmt = 12 ) skip_m(3)
      ELSE IF (ndim == 2) THEN
        skip_m(3) = 2 + 2*nlx + 2*nly
        WRITE( iuni_u, fmt = 11 ) skip_m(3)
      END IF
!
! ... GAS Velocity component X
!
      skip = skip_time + skip_block
      WRITE( iuni_u, fmt = 13 ) 1, filetype, skip
!
! ... SOLID Velocity component X
!
      DO is = 1, nsolid
        skip = skip_time + skip_phase * is + (ngas + 1) * skip_block
        WRITE( iuni_u, fmt = 13 ) is + 1, filetype, skip
      END DO
!
! ... EOF
!
      WRITE( iuni_u, fmt = 14 )
!
!***** V E L O C I T Y   Y ********************************************
! ... Common Header
!
      WRITE( iuni_v, fmt = 1 )
      WRITE( iuni_v, fmt = 2 ) ndim
      IF (ndim == 3) THEN
        WRITE( iuni_v, fmt = 3 ) nx
        WRITE( iuni_v, fmt = 4 ) ny
        WRITE( iuni_v, fmt = 5 ) nz
      ELSE IF (ndim == 2) THEN
        WRITE( iuni_v, fmt = 3 ) nx
        WRITE( iuni_v, fmt = 4 ) nz
      END IF
      WRITE( iuni_v, fmt = 6 ) ndim
      WRITE( iuni_v, fmt = 7 )
      WRITE( iuni_v, fmt = 8 )
!
! ... specify variables and meshes
!
      veclen = nphase
      WRITE( iuni_v, fmt = 9 ) veclen
      skip_m(1) = 2
      WRITE( iuni_v, fmt = 10 ) skip_m(1)
      IF (ndim == 3) THEN
        skip_m(2) = 2 + 2*nlx + nly
        WRITE( iuni_v, fmt = 11 ) skip_m(2)
        skip_m(3) = 2 + 2*nlx + 2*nly
        WRITE( iuni_v, fmt = 12 ) skip_m(3)
      ELSE IF (ndim == 2) THEN
        skip_m(3) = 2 + 2*nlx + 2*nly
        WRITE( iuni_v, fmt = 11 ) skip_m(3)
      END IF
!
! ... GAS Velocity component Y
!
      skip = skip_time + 2 * skip_block
      WRITE( iuni_v, fmt = 13 ) 1, filetype, skip
!
! ... SOLID Velocity component Y
!
      DO is = 1, nsolid
        skip = skip_time + skip_phase * is + (ngas + 2) * skip_block
        WRITE( iuni_v, fmt = 13 ) is + 1, filetype, skip
      END DO
!
! ... EOF
!
      WRITE( iuni_v, fmt = 14 )
!
!***** V E L O C I T Y   Z ********************************************
! ... Common Header
!
      WRITE( iuni_w, fmt = 1 )
      WRITE( iuni_w, fmt = 2 ) ndim
      IF (ndim == 3) THEN
        WRITE( iuni_w, fmt = 3 ) nx
        WRITE( iuni_w, fmt = 4 ) ny
        WRITE( iuni_w, fmt = 5 ) nz
      ELSE IF (ndim == 2) THEN
        WRITE( iuni_w, fmt = 3 ) nx
        WRITE( iuni_w, fmt = 4 ) nz
      END IF
      WRITE( iuni_w, fmt = 6 ) ndim
      WRITE( iuni_w, fmt = 7 )
      WRITE( iuni_w, fmt = 8 )
!
! ... specify variables and meshes
!
      veclen = nphase
      WRITE( iuni_w, fmt = 9 ) veclen
      skip_m(1) = 2 
      WRITE( iuni_w, fmt = 10 ) skip_m(1)
      IF (ndim == 3) THEN
        skip_m(2) = 2 + 2*nlx
        WRITE( iuni_w, fmt = 11 ) skip_m(2)
        skip_m(3) = 2 + 2*nlx + 2*nly + nlz
        WRITE( iuni_w, fmt = 12 ) skip_m(3)
      ELSE IF (ndim == 2) THEN
        skip_m(3) = 2 + 2*nlx + 2*nly + nlz
        WRITE( iuni_w, fmt = 11 ) skip_m(3)
      END IF
!
! ... GAS Velocity component Z
!
      skip = skip_time + ndim * skip_block
      WRITE( iuni_w, fmt = 13 ) 1, filetype, skip
!
! ... SOLID Velocity component Z
!
      DO is = 1, nsolid
        skip = skip_time + skip_phase * is + (ngas + ndim) * skip_block
        WRITE( iuni_w, fmt = 13 ) is + 1, filetype, skip
      END DO
!
! ... EOF
!
      WRITE( iuni_w, fmt = 14 )

!
!**********************************************************************
! ... Common Header
!
  1   FORMAT('# AVS field file')
  2   FORMAT('ndim=',I3)
  3   FORMAT('dim1=',I3)
  4   FORMAT('dim2=',I3)
  5   FORMAT('dim3=',I3)
  6   FORMAT('nspace=',I3)
  7   FORMAT('data=float')
  8   FORMAT('field=rectilinear')

! ... Specify meshes
!
  9   FORMAT('veclen=',I3)
 10   FORMAT('coord 1 file=mesh.dat, filetype=ascii, skip= ',I6)
 11   FORMAT('coord 2 file=mesh.dat, filetype=ascii, skip= ',I6)
 12   FORMAT('coord 3 file=mesh.dat, filetype=ascii, skip= ',I6)

! ... Specify variables
!
 13   FORMAT('variable',I3,' file=./output, filetype = ',A15,' skip= ',I12)
 15   FORMAT('variable',I3,' file=./m_outp, filetype = ',A15,' skip= ',I12)

! ... Common EOF
!
 14   FORMAT('')

      CLOSE( UNIT=iuni_scalar )
      CLOSE( UNIT=iuni_u )
      CLOSE( UNIT=iuni_v )
      CLOSE( UNIT=iuni_w )

      RETURN
      END SUBROUTINE write_AVS_files
!----------------------------------------------------------------------
      SUBROUTINE write_XML_files
      !
      USE control_flags, ONLY: job_type, formatted_output
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions
      USE iotk_module
      USE grid
      USE input_module, ONLY: run_name
      USE time_parameters, ONLY: dt, timestart
      USE kinds
      USE io_files, ONLY: iunxml, xmlfile

      IMPLICIT NONE
 
      CHARACTER(LEN=15) :: filetype
      INTEGER :: skip_m(3)
      INTEGER :: skip, skip_time, skip_block, byte_record
      INTEGER :: lbl, ig, is, it
      INTEGER :: nfields, ndim, veclen, nlx, nly, nlz
      CHARACTER(LEN=80) :: attr
      CHARACTER(LEN=80) :: fldn
      CHARACTER(LEN=11) :: skip_what
      REAL(sgl) :: rsgl

      OPEN( UNIT=iunxml, FILE=xmlfile, STATUS='UNKNOWN')
!
! ... control parameters
!
      IF (job_type == JOB_TYPE_3D) THEN
        ndim = 3
      ELSE IF(job_type == JOB_TYPE_2D) THEN
        ndim = 2
      ELSE
        WRITE(*,*) 'Unknown job_type'
        STOP
      END IF
!
      IF (formatted_output) THEN
        filetype = 'ascii'
      ELSE
        filetype = 'binary'
      END IF

! ... Total number of scalar fields
!
      nfields = nphase*2 + ngas
!
      IF (formatted_output) THEN 
        ! ... number of lines for each block
        lbl = nz * (INT(nx*ny / 10) + MIN(1,MOD(nx*ny,10))) + 4
      ELSE
        ! ... number of bytes (chars) for each block
        lbl = ntot * 4 + 8
      END IF 

! ... Number of lines in mesh file
!
        nlx = ( nx/5 + MIN(1,MOD(nx,5)) + 1 )
        nly = ( ny/5 + MIN(1,MOD(ny,5)) + 1 )
        nlz = ( nz/5 + MIN(1,MOD(nz,5)) + 1 )

      IF (formatted_output) THEN
        skip_time = 4 + 4
        skip_block = lbl
        skip_what = "linestoskip"
      ELSE
        skip_time = 16
        skip_block = lbl
        skip_what = "charstoskip"
      END IF

!***** S C A L A R S **************************************************
! ... Common Header
!
      WRITE( iunxml, fmt="(A39)" ) '<?xml version="1.0" encoding="UTF-8" ?>'

      attr = ' '
      CALL iotk_write_attr( attr, "name", run_name )
      CALL iotk_write_begin( iunxml, "dataset", attr )

        CALL iotk_write_begin( iunxml, "grid" )
          !
          CALL iotk_write_begin( iunxml, "structure" )
            !
            !  Grid dimensions
            ! 
            CALL iotk_write_begin( iunxml, "x" )
              CALL iotk_write_dat( iunxml, "dim", nx )
              skip = 2
              attr = ' '
              CALL iotk_write_attr( attr, "type", "ascii" )
              CALL iotk_write_attr( attr, "linestoskip", skip )
              WRITE( iunxml, * ) '<cfile ' // TRIM(attr) // ' >' // 'mesh.dat' // '</cfile>'
              skip = skip + nlx 
              attr = ' '
              CALL iotk_write_attr( attr, "type", "ascii" )
              CALL iotk_write_attr( attr, "linestoskip", skip )
              WRITE( iunxml, * ) '<file ' // TRIM(attr) // ' >' // 'mesh.dat' // '</file>'
            CALL iotk_write_end( iunxml, "x" )
            !
            IF( ndim == 3 ) THEN
              CALL iotk_write_begin( iunxml, "y" )
                CALL iotk_write_dat( iunxml, "dim", ny )
                skip = skip + nlx 
                attr = ' '
                CALL iotk_write_attr( attr, "type", "ascii" )
                CALL iotk_write_attr( attr, "linestoskip", skip )
                WRITE( iunxml, * ) '<cfile ' // TRIM(attr) // ' >' // 'mesh.dat' // '</cfile>'
                skip = skip + nly
                attr = ' '
                CALL iotk_write_attr( attr, "type", "ascii" )
                CALL iotk_write_attr( attr, "linestoskip", skip )
                WRITE( iunxml, * ) '<file ' // TRIM(attr) // ' >' // 'mesh.dat' // '</file>'
              CALL iotk_write_end( iunxml, "y" )
            ELSE
              skip = skip + nlx
              skip = skip + nly
            END IF
            !
            CALL iotk_write_begin( iunxml, "z" )
              CALL iotk_write_dat( iunxml, "dim", nz )
              skip = skip + nly
              attr = ' '
              CALL iotk_write_attr( attr, "type", "ascii" )
              CALL iotk_write_attr( attr, "linestoskip", skip )
              WRITE( iunxml, * ) '<cfile ' // TRIM(attr) // ' >' // 'mesh.dat' // '</cfile>'
              skip = skip + nlz
              attr = ' '
              CALL iotk_write_attr( attr, "type", "ascii" )
              CALL iotk_write_attr( attr, "linestoskip", skip )
              WRITE( iunxml, * ) '<file ' // TRIM(attr) // ' >' // 'mesh.dat' // '</file>'
            CALL iotk_write_end( iunxml, "z" )
            ! 
          CALL iotk_write_end( iunxml, "structure" )
          !
          CALL iotk_write_begin( iunxml, "attributes" )
            !
            attr = ' '
            CALL iotk_write_attr( attr, "name", "profile" )
            CALL iotk_write_begin( iunxml, "profile", attr )
              !
              attr = ' '
              CALL iotk_write_attr( attr, "type", "ascii" )
              WRITE( iunxml, * ) '<file ' // TRIM(attr) // ' >' // 'improfile.dat' // '</file>'
              !
            CALL iotk_write_end( iunxml, "profile" )
            !
            ! Gas pressure
            !
            skip = skip_time  ! skip the time field and the record marker
            attr = ' '
            CALL iotk_write_attr( attr, "name", "P" )
            CALL iotk_write_begin( iunxml, "scalar", attr )
              attr = ' '
              CALL iotk_write_attr( attr, "type", filetype )
              CALL iotk_write_attr( attr, skip_what, skip )
              WRITE( iunxml, * ) '<file ' // TRIM(attr) // ' >' // 'output' // '</file>'
            CALL iotk_write_end( iunxml, "scalar" )
            !
            ! Gas velocities
            !
            attr = ' '
            CALL iotk_write_attr( attr, "name", "VG" )
            CALL iotk_write_begin( iunxml, "staggered", attr )
              skip = skip + skip_block
              CALL iotk_write_begin( iunxml, "x")
                attr = ' '
                CALL iotk_write_attr( attr, "type", filetype )
                CALL iotk_write_attr( attr, skip_what, skip )
                WRITE( iunxml, * ) '<file ' // TRIM(attr) // ' >' // 'output' // '</file>'
              CALL iotk_write_end( iunxml, "x")
              IF( ndim == 3 ) THEN
                skip = skip + skip_block
                CALL iotk_write_begin( iunxml, "y")
                  attr = ' '
                  CALL iotk_write_attr( attr, "type", filetype )
                  CALL iotk_write_attr( attr, skip_what, skip )
                  WRITE( iunxml, * ) '<file ' // TRIM(attr) // ' >' // 'output' // '</file>'
                CALL iotk_write_end( iunxml, "y")
              END IF
              skip = skip + skip_block
              CALL iotk_write_begin( iunxml, "z")
                attr = ' '
                CALL iotk_write_attr( attr, "type", filetype )
                CALL iotk_write_attr( attr, skip_what, skip )
                WRITE( iunxml, * ) '<file ' // TRIM(attr) // ' >' // 'output' // '</file>'
              CALL iotk_write_end( iunxml, "z")
            CALL iotk_write_end( iunxml, "staggered" )
            !
            ! Gas Temperature
            !
            skip = skip + skip_block
            attr = ' '
            CALL iotk_write_attr( attr, "name", "TG" )
            CALL iotk_write_begin( iunxml, "scalar", attr )
              attr = ' '
              CALL iotk_write_attr( attr, "type", filetype )
              CALL iotk_write_attr( attr, skip_what, skip )
              WRITE( iunxml, * ) '<file ' // TRIM(attr) // ' >' // 'output' // '</file>'
            CALL iotk_write_end( iunxml, "scalar" )
            !
            ! Gas Components
            !
            DO ig = 1, ngas
              skip = skip + skip_block
              attr = ' '
              WRITE( fldn, fmt = "(A3,I1)" ) "XGC", ig
              CALL iotk_write_attr( attr, "name", TRIM( fldn ) )
              CALL iotk_write_begin( iunxml, "scalar", attr )
                attr = ' '
                CALL iotk_write_attr( attr, "type", filetype )
                CALL iotk_write_attr( attr, skip_what, skip )
                WRITE( iunxml, * ) '<file ' // TRIM(attr) // ' >' // 'output' // '</file>'
              CALL iotk_write_end( iunxml, "scalar" )
            END DO
            !
            ! Solid Phases
            !
            DO is = 1, nsolid
              !
              ! Solid Density
              !
              skip = skip + skip_block
              attr = ' '
              WRITE( fldn, fmt = "(A3,I1)" ) "EPS", is
              CALL iotk_write_attr( attr, "name", TRIM( fldn ) )
              CALL iotk_write_begin( iunxml, "scalar", attr )
                attr = ' '
                CALL iotk_write_attr( attr, "type", filetype )
                CALL iotk_write_attr( attr, skip_what, skip )
                WRITE( iunxml, * ) '<file ' // TRIM(attr) // ' >' // 'output' // '</file>'
              CALL iotk_write_end( iunxml, "scalar" )
              !
              ! Solid Velocities
              !
              attr = ' '
              WRITE( fldn, fmt = "(A3,I1)" ) "VS", is
              CALL iotk_write_attr( attr, "name", TRIM( fldn ) )
              CALL iotk_write_begin( iunxml, "staggered", attr )
                skip = skip + skip_block
                CALL iotk_write_begin( iunxml, "x")
                  attr = ' '
                  CALL iotk_write_attr( attr, "type", filetype )
                  CALL iotk_write_attr( attr, skip_what, skip )
                  WRITE( iunxml, * ) '<file ' // TRIM(attr) // ' >' // 'output' // '</file>'
                CALL iotk_write_end( iunxml, "x")
                IF( ndim == 3 ) THEN
                  skip = skip + skip_block
                  CALL iotk_write_begin( iunxml, "y")
                    attr = ' '
                    CALL iotk_write_attr( attr, "type", filetype )
                    CALL iotk_write_attr( attr, skip_what, skip )
                    WRITE( iunxml, * ) '<file ' // TRIM(attr) // ' >' // 'output' // '</file>'
                  CALL iotk_write_end( iunxml, "y")
                END IF
                skip = skip + skip_block
                CALL iotk_write_begin( iunxml, "z")
                  attr = ' '
                  CALL iotk_write_attr( attr, "type", filetype )
                  CALL iotk_write_attr( attr, skip_what, skip )
                  WRITE( iunxml, * ) '<file ' // TRIM(attr) // ' >' // 'output' // '</file>'
                CALL iotk_write_end( iunxml, "z")
              CALL iotk_write_end( iunxml, "staggered" )
              !
              ! Solid Temperature
              !
              skip = skip + skip_block
              attr = ' '
              WRITE( fldn, fmt = "(A2,I1)" ) "TS", is
              CALL iotk_write_attr( attr, "name", TRIM( fldn ) )
              CALL iotk_write_begin( iunxml, "scalar", attr )
                attr = ' '
                CALL iotk_write_attr( attr, "type", filetype )
                CALL iotk_write_attr( attr, skip_what, skip )
                WRITE( iunxml, * ) '<file ' // TRIM(attr) // ' >' // 'output' // '</file>'
              CALL iotk_write_end( iunxml, "scalar" )
              !
            END DO
            !
          CALL iotk_write_end( iunxml, "attributes" )
          !
        CALL iotk_write_end( iunxml, "grid" )
        !
        attr = ' '
        rsgl = timestart
        CALL iotk_write_attr( attr, "timestart", rsgl )
        rsgl = dt
        CALL iotk_write_attr( attr, "timescale", rsgl )
        CALL iotk_write_attr( attr, "start", first_out )
        CALL iotk_write_attr( attr, "stride", incr_out )
        CALL iotk_write_begin( iunxml, "frames", attr )
!        DO it = first_out, last_out, incr_out
!          CALL iotk_write_dat( iunxml, "timestep", it )
!        END DO
          WRITE( iunxml, * ) (last_out - first_out + 1)/incr_out
        CALL iotk_write_end( iunxml, "frames" )
        !
      CALL iotk_write_end( iunxml, "dataset" )

      CLOSE( UNIT=iunxml )
!
      RETURN
      END SUBROUTINE write_XML_files
!-----------------------------------------------------------------------
      END MODULE postp_output
!----------------------------------------------------------------------
