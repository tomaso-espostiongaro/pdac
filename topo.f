!----------------------------------------------------------------------
      MODULE volcano_topography
!----------------------------------------------------------------------
      USE dimensions, ONLY: nx, ny, nz, ntot
      USE parallel, ONLY: mpime, root

      IMPLICIT NONE
!
! ... This is a standard ascii dem file: the header contains
! ... information about the resolution and the geo-referencing.
! ... The UTM map number is lost.
! ... 'ncols' and 'nrows' are the number of pixel in x and y directions
! ... 'xllcorner' is the distance in metres from the reference meridian
! ...             of the lowerleft corner + 500,000m !
! ... 'yllcorner' is the distance in metres from the equator of the
! ...             lowerleft corner
! ... 'cellsize'  is the cell spacing in metres
! ... 'nodata_value' is the value assigned to missing data 
! ... 'elev'      is the elevation (in centimetres)
!
      TYPE dem_ascii
        INTEGER :: ncols      
        INTEGER :: nrows      
        REAL*8  :: xllcorner  
        REAL*8  :: yllcorner  
        REAL*8  :: cellsize   
        REAL*8  :: nodata_value
        REAL*8, ALLOCATABLE :: elev(:)
      END TYPE dem_ascii
!
! ... coordinates of the input mesh
!
      REAL*8, ALLOCATABLE, DIMENSION(:) :: cx
      REAL*8, ALLOCATABLE, DIMENSION(:) :: cy
      REAL*8, ALLOCATABLE, DIMENSION(:) :: cz
!
! ... arrays used for the interpolation
!
      INTEGER, ALLOCATABLE :: next(:)  ! index of the topographic point
                                       ! located right of a grid point 
      INTEGER, ALLOCATABLE :: ord(:)   ! z-coord of grid locations laying
                                       ! below the profile
!
! ... dist defines implicitly the profile
!
      REAL*8, ALLOCATABLE  :: dist(:)  ! vertical distance above topography

! ... input topography points
!
      REAL*8, ALLOCATABLE :: xtop(:), ytop(:), ztop(:) 
      INTEGER :: noditop
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE import_topo
      USE control_flags, ONLY: job_type
!
      IMPLICIT NONE
!
! ... Read the topography file
!
      IF (job_type == '2D') THEN
        CALL read_2Dprofile
      ELSE IF (job_type == '3D') THEN
        CALL import_dem_ascii
        !CALL compute_UTM_coord
      END IF
!
      ALLOCATE(next(nx))
      ALLOCATE(ord(nx))
      ALLOCATE(dist(ntot))
!
      DEALLOCATE (xtop, ztop)
      DEALLOCATE (cx, cy, cz)
      DEALLOCATE (next, ord, dist)
!
      RETURN
      END SUBROUTINE import_topo
!----------------------------------------------------------------------
      SUBROUTINE read_2Dprofile
      USE grid, ONLY: topography
      IMPLICIT NONE

      INTEGER :: n

      REAL*8 :: ratio
      REAL*8 :: xmin, xmax, zmin, zmax
      CHARACTER(LEN=20) :: topo_file
!
! ... read the topography, as defined on a generic rectilinear
! ... mesh "xtop"
!
      topo_file = TRIM(topography)
      OPEN(UNIT=3, FILE=topo_file, STATUS='OLD')

      READ(3,*) noditop

      ALLOCATE(xtop(noditop))
      ALLOCATE(ztop(noditop))

      DO n=1, noditop
        READ(3,*) xtop(n),ztop(n)
      END DO

      CLOSE(3)
!
      xmin = MINVAL(xtop)
      xmax = MAXVAL(xtop)
      zmin = MINVAL(ztop)
      zmax = MAXVAL(ztop)

      ratio=(xmax-xmin)/(zmax-zmin)

      END SUBROUTINE read_2Dprofile
!----------------------------------------------------------------------
      SUBROUTINE import_dem_ascii
      USE grid, ONLY: topography
      IMPLICIT NONE

      CHARACTER(LEN=20) :: topo_file
      TYPE(dem_ascii) :: vdem
      INTEGER :: nc, nr, p

      topo_file = TRIM(topography)
      OPEN(UNIT=3, FILE=topo_file, STATUS='OLD')

      READ(3,*) vdem%ncols
      READ(3,*) vdem%nrows
      READ(3,*) vdem%xllcorner
      READ(3,*) vdem%yllcorner
      READ(3,*) vdem%cellsize
      READ(3,*) vdem%nodata_value

      nc = vdem%ncols
      nr = vdem%nrows
      ALLOCATE( vdem%elev(nc*nr) )

      READ(3,*) ( vdem%elev(p), p=1,nc*nr )
!
      RETURN
      END SUBROUTINE import_dem_ascii
!----------------------------------------------------------------------
!      SUBROUTINE compute_UTM_coordinates
!      IMPLICIT NONE
!      
!      RETURN
!      END SUBROUTINE compute_UTM_coordinates
!----------------------------------------------------------------------
      SUBROUTINE interpolate_2d(topo,ff)
!
! ... interpolate the topographic profile on a given
! ... computational mesh (either centered or staggered)

      IMPLICIT NONE

      LOGICAL, INTENT(OUT), DIMENSION(:) :: ff
      REAL*8, INTENT(OUT), DIMENSION(:) :: topo

      INTEGER :: i, k, ijk, l, n
      REAL*8 :: grad
!
! ... locate the grid points near the topographic points
! ... and interpolate linearly the profile  
!
      dist = 1.0D10
      ff = .FALSE.
!
! ... locate the topography on the mesh
!
      l=1
      DO n = 1, noditop
        DO i = l, nx

          IF (xtop(n) >= cx(i)) THEN

            ! ... next(i) indicates the progressive number of the
            ! ... topographic point laying on the right of a grid center 'i'
            ! ... 'l' counts the grid points
            !
            next(i) = n
 
            IF (n == 1) THEN
              topo(i) = ztop(1)
            ELSE
              grad = (ztop(n)-ztop(n-1))/(xtop(n)-xtop(n-1))
              topo(i) = ztop(n-1) + (cx(i)-xtop(n-1)) * grad
            END IF

	    l=l+1

            DO k = 1, nz
              ! ... dist defines implicitly the profile
              !
	      ijk = i + (k-1) * nx
	      dist(ijk) = cz(k) - topo(i)
	    ENDDO

          ENDIF
        ENDDO
      ENDDO
!
! ... ord is the last mesh point laying below the topography
!
      DO i = 1, nx
        DO k = 1, nz
          IF (cz(k) <= topo(i))  ord(i) = k
        END DO
      END DO
!
! ... locate forced point on the mesh (more than one forcing points can
! ... lay on the same coordinate on the plane)
!
      ijk = 1 + (ord(1)-1)*nx
      ff(ijk) = .TRUE.
      DO i = 2, nx
        IF (ord(i) >= ord(i-1)) THEN
          DO k = ord(i-1), ord(i)
            ijk = i + (k-1) * nx
            ff(ijk) = .TRUE.
          END DO
        ELSE IF ( ord(i) < ord(i-1) ) THEN
          ijk = i + (ord(i)-1) * nx
          ff(ijk) = .TRUE.
          DO k = ord(i), ord(i-1)
            ijk = (i-1) + (k-1) * nx
            ff(ijk) = .TRUE.
          END DO
        END IF
      END DO
!
      RETURN
      END SUBROUTINE interpolate_2d
!----------------------------------------------------------------------
      SUBROUTINE grid_locations(sx,sy,sz)
      USE grid, ONLY: x, xb, y, yb, z, zb
      USE dimensions, ONLY: nx, ny, nz

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: sx, sy, sz
!
! ... initialize the grids
!
      IF(.NOT.(ALLOCATED(cx))) ALLOCATE(cx(nx))
      IF(.NOT.(ALLOCATED(cy))) ALLOCATE(cy(ny))
      IF(.NOT.(ALLOCATED(cz))) ALLOCATE(cz(nz))
      cx = 0.D0
      cy = 0.D0
      cz = 0.D0
!
      IF( sx == 1 ) THEN

        cx = xb
        cy = y
        cz = z

      ELSE IF( sy == 1 ) THEN

        cx = x
        cy = yb
        cz = z

      ELSE IF( sz == 1 ) THEN

        cx = x
        cy = y
        cz = zb
        
      ELSE 
        
        cx = x
        cy = y
        cz = z

      END IF

      END SUBROUTINE grid_locations
!----------------------------------------------------------------------
      END MODULE volcano_topography
!----------------------------------------------------------------------
