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
! ... 'nx' and 'ny' are the number of pixel in x and y directions
! ... 'xcorner' is the distance in metres from the reference meridian
! ...             of the upperleft corner + 500,000m !
! ... 'ycorner' is the distance in metres from the equator of the
! ...             upperleft corner
! ... 'cellsize'  is the cell spacing in metres
! ... 'nodata_value' is the value assigned to missing data 
! ... 'elev'      is the elevation (in centimetres)
!
      TYPE dem_header
        INTEGER :: nx      
        INTEGER :: ny      
        REAL*8  :: xcorner  
        REAL*8  :: ycorner  
        REAL*8  :: cellsize   
        REAL*8  :: nodata_value
      END TYPE dem_header

      TYPE(dem_header) :: vdem
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
      INTEGER, ALLOCATABLE :: nextx(:)  ! index of the topographic point
                                        ! located right of a grid point 
      INTEGER, ALLOCATABLE :: nexty(:)  ! index of the topographic point
                                        ! located right of a grid point 
      INTEGER, ALLOCATABLE :: ord(:)   ! z-coord of grid locations laying
                                       ! below the profile
      INTEGER, ALLOCATABLE :: ord2d(:,:)   ! z-coord of grid locations laying
                                           ! below the profile
!
! ... 'dist' defines implicitly the profile
!
      REAL*8, ALLOCATABLE  :: dist(:)  ! vertical distance above topography

! ... input topography points
!
      REAL*8, ALLOCATABLE :: xtop(:), ytop(:), ztop(:), ztop2d(:,:)
      INTEGER :: noditop, noditopx, noditopy
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE import_topography
!
! ... Read the topography file
!
      USE control_flags, ONLY: job_type
      IMPLICIT NONE

      IF (job_type == '2D') THEN
        CALL read_2Dprofile
      ELSE IF (job_type == '3D') THEN
        CALL read_dem_ascii
        CALL compute_UTM_coords
      END IF

      CALL set_profile
!
      RETURN
      END SUBROUTINE import_topography
!----------------------------------------------------------------------
      SUBROUTINE read_2Dprofile
      USE grid, ONLY: topography
      USE parallel, ONLY: mpime, root
      IMPLICIT NONE

      INTEGER :: n

      REAL*8 :: ratio
      REAL*8 :: xmin, xmax, zmin, zmax
      CHARACTER(LEN=80) :: topo_file
!
! ... read the topography, as defined on a generic rectilinear
! ... mesh "xtop"
!
      topo_file = TRIM(topography)

      IF (mpime == root) THEN
        OPEN(UNIT=3, FILE=topo_file, STATUS='OLD')
        READ(3,*) noditop
      END IF

      CALL bcast_integer(noditop,1,root)

      ALLOCATE(xtop(noditop))
      ALLOCATE(ztop(noditop))

      IF (mpime == root) THEN
        DO n=1, noditop
          READ(3,*) xtop(n),ztop(n)
        END DO
        CLOSE(3)
      END IF

      CALL bcast_real(xtop,noditop,root)
      CALL bcast_real(ytop,noditop,root)
!
      xmin = MINVAL(xtop)
      xmax = MAXVAL(xtop)
      zmin = MINVAL(ztop)
      zmax = MAXVAL(ztop)

      ratio=(xmax-xmin)/(zmax-zmin)

      END SUBROUTINE read_2Dprofile
!----------------------------------------------------------------------
      SUBROUTINE read_dem_ascii
!
! ... Read the standard ASCII DEM format
!
      USE grid, ONLY: topography
      USE parallel, ONLY: mpime, root
      IMPLICIT NONE

      CHARACTER(LEN=80) :: topo_file
      INTEGER :: p
      INTEGER :: i, j
      
      REAL    :: xll, yll
      REAL    :: dd
      INTEGER :: noval
      INTEGER :: elevation

      IF (mpime == root) THEN
        topo_file = TRIM(topography)
        OPEN(UNIT=3, FILE=topo_file, STATUS='OLD')
  
        WRITE(*,*) 'Reading topography file: ', topo_file

        READ(3,*) noditopx
        READ(3,*) noditopy
        READ(3,*) xll
        READ(3,*) yll
        READ(3,*) dd
        READ(3,*) noval

      END IF

      CALL bcast_integer(noditopx,1,root)
      CALL bcast_integer(noditopy,1,root)
      CALL bcast_real(xll,1,root)
      CALL bcast_real(yll,1,root)
      CALL bcast_real(dd,1,root)
      CALL bcast_integer(noval,1,root)

      vdem%nx           = noditopx
      vdem%ny           = noditopy
      vdem%xcorner      = xll
      vdem%ycorner      = yll
      vdem%cellsize     = dd
      vdem%nodata_value = noval
!
      ALLOCATE(ztop2d(vdem%nx,vdem%ny))
!
      ALLOCATE(xtop(vdem%nx))
      ALLOCATE(ytop(vdem%ny))
!
      DO i = 1, vdem%nx
        xtop(i) = vdem%xcorner + (i-1) * vdem%cellsize
      END DO

      DO j = vdem%ny, 1, -1
        ytop(j) = vdem%ycorner - (vdem%ny - j) * vdem%cellsize
      END DO

      IF (mpime == root) THEN
        DO j = vdem%ny, 1, -1
          DO i = 1, vdem%nx
            READ(3,*) elevation
            ztop2d(i,j) = DBLE(elevation) / 100.D0
          END DO
        END DO
      END IF

      CALL bcast_real(ztop2d,noditopx*noditopy,root)
!
      RETURN
      END SUBROUTINE read_dem_ascii
!----------------------------------------------------------------------
      SUBROUTINE compute_UTM_coords
      USE grid, ONLY: center_x, center_y
      USE grid, ONLY: x, y, xb, yb
      USE grid, ONLY: iv, jv
      IMPLICIT NONE

      REAL*8 :: transl_x, transl_y

      transl_x = center_x - x(iv)
      transl_y = center_y - y(jv)

      x  = x  + transl_x
      xb = xb + transl_x
      y  = y  + transl_y
      yb = yb + transl_y
!
      OPEN(17,FILE='mesh.dat')
!
      WRITE(17,*) 'Geoferenced x-y mesh'

      WRITE(17,*) 'x'
      WRITE(17,*) x
      WRITE(17,*) 'xb'
      WRITE(17,*) xb
      WRITE(17,*) 'y'
      WRITE(17,*) y
      WRITE(17,*) 'yb'
      WRITE(17,*) yb
      
      WRITE(6,*) 'mesh center: ', iv, jv
      WRITE(6,*) 'center coordinates: ', x(iv), y(jv)
      
      RETURN
      END SUBROUTINE compute_UTM_coords
!----------------------------------------------------------------------
      SUBROUTINE vertical_shift(topo2d)
      USE grid, ONLY: z, zb, dz, dzmin
      USE grid, ONLY: iv, jv, kv
      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:,:) :: topo2d
      REAL*8 :: transl_z

      !
      ! ... Translate vertically the numerical mesh to minimize the
      ! ... number of topographic cells
      !
      transl_z = MINVAL(topo2d)
      WRITE(*,*) 'Minimum topographic quota: ', transl_z

      z  = z  + transl_z
      zb = zb + transl_z

      WRITE(17,*) 'z'
      WRITE(17,*) z
      WRITE(17,*) 'zb'
      WRITE(17,*) zb

      CLOSE(17)
!
! ... define the ordinate of the vent
!
      kv = ord2d(iv,jv)

      RETURN
      END SUBROUTINE vertical_shift
!----------------------------------------------------------------------
      SUBROUTINE interpolate_2d(topo,ff)
!
! ... interpolate the topographic profile on a given
! ... computational mesh (either centered or staggered)
! ... as defined by the arrays (cx, cy, cz)

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
! ... interpolate the topography on the mesh
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
      SUBROUTINE interpolate_dem(topo2d, ff)
      USE dimensions, ONLY: nx, ny, nz
      IMPLICIT NONE

      REAL*8, INTENT(OUT), DIMENSION(:,:) :: topo2d
      LOGICAL, INTENT(OUT), DIMENSION(:) :: ff

      REAL*8 xmin,xmax,ymin,ymax,zmin,zmax !range dei valori della top.
      REAL*8 ratio

      REAL*8 dist1y,dist2y,dist1x,dist2x,alpha,beta
      INTEGER i,j,k,h,l,ii,jj
      INTEGER ijk, tp1, tp2

      ff = .FALSE.
        
!C============================================
!C===    trova le posizioni dei nodi      ====
!C===    della nuova griglia rispetto     ====
!C===    alla griglia iniziale            ====
!C============================================

      l=1
      DO i = 1, noditopx
        DO ii = l, nx
      
!============================================
!===    cerca i nodi della griglia che    ===
!===    che stanno a sx. di xtop(i)       ===
!============================================
    
          IF (xtop(i) >= cx(ii)) THEN
            nextx(ii)=i
            l=l+1
          ENDIF
        ENDDO
      ENDDO
	
      l=1
      DO j = 1, noditopy
        DO jj = l, ny

!C============================================
!C===    cerca i nodi della griglia che    ===
!C===    che stanno a sotto  ytop(i)       ===
!C============================================

          IF (ytop(j) >= cy(jj)) THEN
            nexty(jj)=j
            l=l+1
          ENDIF
        ENDDO
      ENDDO


! il nodo della nuova griglia di indici (i,j) sara' allora
! contenuto nel rettangolo con i vertici con indici:
! P1=nextx(i),nexty(j)
! P2=nextx(i-1),nexty(j)
! P3=nextx(i-1),nexty(j-1)
! P4=nextx(i),nexty(j-1)
	

! sulla nuova griglia interpoliamo le quote di input ztop per 
! ottenere la quota coorZ nel punto di indici (i,j)
 

! interpolazione bilineare sui nodi interni (1<i<nodiGRIDx, 1<j<nodigGRIDy)
! utilizzando le quote nei punti P1,..,P4 definiti sopra

	DO i=1,nx

	   dist1x = cx(i) - xtop(nextx(i)-1)
	   dist2x = xtop(nextx(i)) - cx(i)
	   alpha  = dist1x/(dist1x+dist2x)

	   DO j=1,ny
           
	      dist1y = cy(j) - ytop(nexty(j)-1)
	      dist2y = ytop(nexty(j)) - cy(j)
	      beta   = dist1y/(dist1y+dist2y)

              tp1    = alpha * ztop2d(nextx(i),nexty(j))   + &
                       (1.D0 - alpha) * ztop2d(nextx(i)-1,nexty(j))

              tp2    = alpha * ztop2d(nextx(i),nexty(j)-1) + &
                       (1.D0 - alpha) * ztop2d(nextx(i)-1,nexty(j)-1)

              topo2d(i,j) = beta * tp1 + (1.D0 - beta) * tp2

	   ENDDO

	ENDDO

! ... Locate cells laying below the topography 'ord2d(i,j)'

	DO j=1,ny
	   DO i=1,nx

	      DO k=1,nz

		IF (cz(k) <= topo2d(i,j)) THEN
		   ord2d(i,j) = k  
		ENDIF

	      ENDDO

	   ENDDO
	ENDDO
!
! ... Identify forcing points
! ... (check boundaries!)
!
      i = 1
      DO j = 1, ny
        k = ord2d(i,j)
        ijk = i + (j-1) * nx + (k-1) * nx * ny
        ff(ijk) = .TRUE.
      END DO

      i = nx
      DO j = 1, ny 
        k = ord2d(i,j)
        ijk = i + (j-1) * nx + (k-1) * nx * ny
        ff(ijk) = .TRUE.
      END DO

      DO i = 2, nx - 1

        j = 1
        k = ord2d(i,j)
        ijk = i + (j-1) * nx + (k-1) * nx * ny
        ff(ijk) = .TRUE.

        j = ny
        k = ord2d(i,j)
        ijk = i + (j-1) * nx + (k-1) * nx * ny
        ff(ijk) = .TRUE.

        DO j = 2, ny - 1
          DO k = 1, ord2d(i,j) - 1
            ijk = i + (j-1) * nx + (k-1) * nx * ny
            IF( ( ord2d(i-1,j)   < k ) .OR. &
                ( ord2d(i-1,j-1) < k ) .OR. &
                ( ord2d(i,j-1)   < k ) .OR. &
                ( ord2d(i+1,j-1) < k ) .OR. &
                ( ord2d(i+1,j)   < k ) .OR. &
                ( ord2d(i+1,j+1) < k ) .OR. &
                ( ord2d(i,j+1)   < k ) .OR. &
                ( ord2d(i-1,j+1) < k ) ) ff(ijk) = .TRUE.
          END DO
          
          k = ord2d(i,j)
          ijk = i + (j-1) * nx + (k-1) * nx * ny
          ff(ijk) = .TRUE.
        
        END DO
      END DO
!      
      RETURN
      END SUBROUTINE interpolate_dem
!----------------------------------------------------------------------
      SUBROUTINE grid_locations(sx,sy,sz)
      USE grid, ONLY: x, xb, y, yb, z, zb
      USE dimensions, ONLY: nx, ny, nz

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: sx, sy, sz
!
! ... initialize the grids
!
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
      SUBROUTINE set_profile
!
! ... Import the topography from standard ascii formats.
! ... Interpolate topography on cell centers and write the
! ... implicit profile
!
      USE control_flags, ONLY: job_type, lpr, immb, itp
      USE grid, ONLY: z

      IMPLICIT NONE
      INTEGER :: i,j,k,ijk
      LOGICAL, ALLOCATABLE :: dummy(:)
      REAL*8, ALLOCATABLE  :: topo(:)
      REAL*8, ALLOCATABLE  :: topo2d(:,:)

      ALLOCATE (dummy(ntot))
      ALLOCATE (dist(ntot))
!
      ALLOCATE (cx(nx))
      ALLOCATE (cy(ny))
      ALLOCATE (cz(nz))
!
      IF (job_type == '2D') THEN
  
        ALLOCATE(topo(nx))

        ALLOCATE(next(nx))
        ALLOCATE(ord(nx))

        CALL grid_locations(0,0,0)
        CALL interpolate_2d(topo, dummy)

      ELSE IF (job_type == '3D') THEN

        ALLOCATE(topo2d(nx,ny))

        ALLOCATE(nextx(nx))
        ALLOCATE(nexty(ny))
        ALLOCATE(ord2d(nx,ny))

        CALL grid_locations(0,0,0)
        CALL interpolate_dem(topo2d, dummy)
        CALL vertical_shift(topo2d)

      END IF

      IF (job_type == '2D') THEN
        DEALLOCATE(topo)
      ELSE IF (job_type == '3D') THEN
        DEALLOCATE(topo2d)
      END IF
      DEALLOCATE (dummy)

      RETURN
      END SUBROUTINE set_profile
!----------------------------------------------------------------------
      SUBROUTINE write_profile
      USE control_flags, ONLY: job_type
!
      IMPLICIT NONE

      INTEGER :: i, j, k, ijk

      IF (mpime == root) THEN
        OPEN(UNIT=14,FILE='improfile.dat',STATUS='UNKNOWN')
        
        DO k=1,nz
          DO j=1,ny
            DO i=1,nx
              ! ... dist defines implicitly the profile
              !
              ijk = i + (j-1) * nx + (k-1) * nx * ny
              WRITE(14,*) dist(ijk)
            END DO
          END DO
        END DO

        CLOSE(14)
      END IF
!
! ... set flags on topography
!
      CALL set_flag3
!
! ... Deallocate all arrays in the topography module
!
      IF (job_type == '2D') THEN
        DEALLOCATE (next)
        DEALLOCATE (ord)
        DEALLOCATE (xtop, ztop)
      ELSE IF (job_type == '3D') THEN
        DEALLOCATE(nextx, nexty)
        DEALLOCATE(ord2d)
        DEALLOCATE (xtop, ytop, ztop2d)
      END IF
      DEALLOCATE (dist)
      DEALLOCATE (cx, cy, cz)
! 
      RETURN
      END SUBROUTINE write_profile
!----------------------------------------------------------------------
      SUBROUTINE set_flag3
!
! ... Set cell-flag = 3 in those cells belonging to the topography
! ... where forcing is not applied. Values of fields in these cells
! ... are set to zero when initialized or kept undefined
!
      USE control_flags, ONLY: lpr, job_type
      USE grid, ONLY: fl, z
      IMPLICIT NONE

      INTEGER :: i, j, k, ijk
!
      IF( job_type == '2D') THEN
        DO i=2, nx-1
          DO k = 1, nz
            ijk = i + (k-1) * nx
            IF (ord(i) >= k) THEN
              fl(ijk) = 3
            ELSE
              EXIT
            END IF
          END DO
        END DO
      ELSE IF( job_type == '3D') THEN
        DO j=2, ny-1
          DO i=2, nx-1
            DO k = 1, nz
              ijk = i + (j-1) * nx + (k-1) * nx * ny
              IF (ord2d(i,j) >= k) THEN
                fl(ijk) = 3
              ELSE
                EXIT
              END IF
            END DO
          END DO
        END DO
      END IF
!
      END SUBROUTINE set_flag3
!----------------------------------------------------------------------
      END MODULE volcano_topography
!----------------------------------------------------------------------
