!----------------------------------------------------------------------
      MODULE io_serial
!----------------------------------------------------------------------
      USE dimensions
      USE kinds
!      USE control_flags, ONLY: job_type
!      USE io_files, ONLY: logunit
!      USE output_dump, ONLY: formatted_output
!
      IMPLICIT NONE
!
! ... main fields
!
      REAL*8, ALLOCATABLE :: improfile(:,:,:)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: array_map
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: topo2d

      INTEGER :: first_out, last_out, incr_out
      REAL*8 :: deltaz
! ... not used
!      INTEGER :: variable_n, field_n

      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE read_array( iunit, sarray, lform )

      !  This subroutine reads a REAL array ( sarray ) of
      !  ntot elements from file iunit
      !  iunit  (input)  file to be read
      !  array  (output) data read from file and distributed to processors
      !  lform  (input)  format of the file 
      !                  .TRUE.  = formatted
      !                  .FALSE. = unformatted

      INTEGER, INTENT(IN) :: iunit
      LOGICAL, INTENT(IN) :: lform
      REAL :: sarray(:)

      INTEGER :: ijk, ierr

      IF( lform ) THEN
            READ(iunit,*) ( sarray(ijk), ijk = 1, ntot )
      ELSE
            READ(iunit) ( sarray(ijk), ijk = 1, ntot )
      END IF

      RETURN
      END SUBROUTINE read_array
!----------------------------------------------------------------------
      SUBROUTINE write_array( iunit, sarray, lform )

      !  This subroutine writes a REAL array ( sarray ) of
      !  ntot elements from file iunit
      !  iunit  (input)  file to be read
      !  array  (output) data read from file and distributed to processors
      !  lform  (input)  format of the file 
      !                  .TRUE.  = formatted
      !                  .FALSE. = unformatted

      INTEGER, INTENT(IN) :: iunit
      LOGICAL, INTENT(IN) :: lform
      REAL :: sarray(:)

      INTEGER :: ijk, ierr

      IF( lform ) THEN
        WRITE(iunit,100) ( sarray(ijk), ijk = 1, ntot )
      ELSE
        WRITE(iunit) ( sarray(ijk), ijk = 1, ntot )
      END IF

 100  FORMAT( 5(G14.6E3,1X) )
      RETURN
      END SUBROUTINE write_array
!----------------------------------------------------------------------
      SUBROUTINE read_implicit_profile
      USE grid, ONLY: x, y, z, xb, yb, zb
      USE io_files, ONLY: tempunit

      !  This subroutine reads the implicit profile
      !  and the mesh file computed by PDAC when a 3D
      !  volcano topography is imported from a DEM file

      INTEGER :: i,j,k
      ALLOCATE(improfile(nx,ny,nz))
!
! ... Read the georeferenced mesh
!
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
!      
! ... Read the Implicit Profile
!
      OPEN(tempunit,FILE='improfile.dat',STATUS='OLD')
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            READ(tempunit,*) improfile(i,j,k)
          END DO
        END DO
      END DO
      CLOSE(tempunit)

 100  FORMAT(5(F20.6))
      RETURN
      END SUBROUTINE read_implicit_profile
!-----------------------------------------------------------------------
      SUBROUTINE write_topo2d

      USE control_flags, ONLY: job_type
      USE dimensions, ONLY: nx, ny, nz
      USE grid, ONLY: z
      USE io_files, ONLY: tempunit

      IMPLICIT NONE

      REAL*8 :: quota
      INTEGER :: i, j, k, ijk

      IF (job_type == '2D') RETURN

      ALLOCATE(topo2d(nx,ny))
      topo2d = -9999 

      DO i = 1, nx
        DO j = 1, ny
          !
          search: DO k = 1, nz
            quota = improfile(i,j,k)
            IF (quota >= 0.D0 .AND. topo2d(i,j) == -9999) THEN
              topo2d(i,j) = z(k) - quota
              EXIT search
            END IF
          END DO search
          !
        END DO
      END DO
!
! ... Print out the new 2D DEM file
!
      OPEN(UNIT=tempunit,FILE='topo2d.dat')
      DO j = 1, ny
          WRITE(tempunit,122) (topo2d(i,j), i=1, nx)
      END DO
      CLOSE(tempunit)

 122  FORMAT(10(1x,G14.6E3))

      RETURN
      END SUBROUTINE write_topo2d
!-----------------------------------------------------------------------
      SUBROUTINE write_map(nfil,array,labl)

      USE control_flags, ONLY: job_type
      USE dimensions, ONLY: nx, ny, nz
      USE grid, ONLY: z
      USE io_files, ONLY: tempunit
      IMPLICIT NONE

      REAL, INTENT(IN), DIMENSION(:) :: array
      INTEGER, INTENT(IN) :: nfil
      CHARACTER(LEN=2), INTENT(IN) :: labl
      REAL*8 :: alpha, map, quota
      INTEGER :: i, j, k, ijk, ijkm
      CHARACTER( LEN = 4 ) :: lettera
      CHARACTER( LEN = 20 ) :: filnam

      IF (job_type == '2D') RETURN

      ALLOCATE(array_map(nx,ny))
      array_map = 0.D0
      
      filnam='map_'//labl//'.'//lettera(nfil)
!
      alpha = 0.0D0
      DO i = 1, nx
        DO j = 1, ny
          !
          search: DO k = 1, nz
            quota = improfile(i,j,k)
            IF (quota >= deltaz) THEN
                ijk  = i + (j-1) * nx + (k-1) * nx * ny
                ijkm = i + (j-1) * nx + (k-2) * nx * ny
                alpha = deltaz - improfile(i,j,k-1)
                alpha = alpha / (z(k) - z(k-1))
                !
                !...errore!!!
                !map = alpha* array(ijk) + (1.D0-alpha) * array(ijkm)
                map = alpha* array(ijk) + (1.D0-alpha) * array(ijk)
                ! ... Map the value reached at any given position at
                ! ... given time
                array_map(i,j) = map
                EXIT search
            END IF
          END DO search
          !
        END DO
      END DO
!
! ... Print out the map and the new 2D DEM file
!
      OPEN(UNIT=tempunit,FILE=filnam)
      DO j = 1, ny
          WRITE(tempunit,122) (array_map(i,j), i=1, nx)
      END DO
      CLOSE(tempunit)

 122  FORMAT(10(1x,G14.6E3))

      DEALLOCATE(array_map)

      RETURN
      END SUBROUTINE write_map
!----------------------------------------------------------------------
      END MODULE io_serial
!----------------------------------------------------------------------
