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
! ... 'xcorner' is the distance in metres from the reference meridian
! ...             of the lowerleft corner + 500,000m !
! ... 'ycorner' is the distance in metres from the equator of the
! ...             lowerleft corner
! ... 'cellsize'  is the cell spacing in metres
! ... 'nodata_value' is the value assigned to missing data 
! ... 'elev'      is the elevation (in centimetres)
!
      TYPE dem_ascii
        INTEGER :: ncols      
        INTEGER :: nrows      
        REAL*8  :: xcorner  
        REAL*8  :: ycorner  
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
      REAL*8, ALLOCATABLE :: xtop(:), ytop(:), ztop(:) 
      REAL*8, ALLOCATABLE :: ztop2d(:,:) 
      INTEGER :: noditop, noditopx, noditopy
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE read_topo
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
!
      RETURN
      END SUBROUTINE read_topo
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
      SUBROUTINE read_dem_ascii
!
! ... Read the standard ASCII DEM format
!
      USE grid, ONLY: topography
      IMPLICIT NONE

      CHARACTER(LEN=20) :: topo_file
      TYPE(dem_ascii) :: vdem
      INTEGER :: p
      INTEGER :: i, j

      topo_file = TRIM(topography)
      OPEN(UNIT=3, FILE=topo_file, STATUS='OLD')

      READ(3,*) vdem%ncols
      READ(3,*) vdem%nrows
      READ(3,*) vdem%xcorner
      READ(3,*) vdem%ycorner
      READ(3,*) vdem%cellsize
      READ(3,*) vdem%nodata_value

      noditopx = vdem%ncols
      noditopy = vdem%nrows
      ALLOCATE( vdem%elev(noditopx*noditopy) )

      DO p = 1,noditopx*noditopy
        READ(3,*) vdem%elev(p)
        vdem%elev(p) = vdem%elev(p) / 100.D0 ! elevation in metres
      END DO
!
      ALLOCATE(xtop(noditopx))
      ALLOCATE(ytop(noditopy))
      ALLOCATE(ztop2d(noditopx,noditopy))
!
      DO i = 1, noditopx
        xtop(i) = vdem%xcorner + (i-1) * vdem%cellsize
      END DO

      DO j = noditopy, 1, -1
        ytop(j) = vdem%ycorner - (noditopy - j) * vdem%cellsize

        DO i = 1, noditopx
          ztop2d(i,j) = vdem%elev( i + (noditopy - j) * noditopx )
        END DO

      END DO
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

!      !!!WARNING!!! controlla traslazione a.s.l. !!!
!      WRITE(*,*) 'Minimum topographic quota: ', MINVAL(ztop2d)
!      transl_z = MINVAL(ztop2d) - z(kv) 
!      z  = z  + transl_z
!      zb = zb + transl_z
!      WRITE(17,*) z
!      WRITE(17,*) zb

      x  = x  + transl_x
      xb = xb + transl_x
      y  = y  + transl_y
      yb = yb + transl_y

      WRITE(17,*) 'Geoferenced x-y mesh'

      WRITE(17,*) x
      WRITE(17,*)
      WRITE(17,*) y
      WRITE(17,*)
      WRITE(17,*) xb
      WRITE(17,*)
      WRITE(17,*) yb
      
      WRITE(6,*) 'mesh center: ', iv, jv
      WRITE(6,*) 'center coordinates: ', x(iv), y(jv)
      
      RETURN
      END SUBROUTINE compute_UTM_coords
!----------------------------------------------------------------------
      SUBROUTINE interpolate_2d(topo,ff)
!
! ... interpolate the topographic profile on a given
! ... computational mesh (either centered or staggered)
! ... as defined by the arrays (cx, cy, cz)

      IMPLICIT NONE

      LOGICAL, INTENT(OUT), DIMENSION(:), OPTIONAL :: ff
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
      SUBROUTINE interpolate_dem(topo2d)
      USE dimensions, ONLY: nx, ny, nz
      IMPLICIT NONE
      REAL*8, INTENT(OUT), DIMENSION(:,:) :: topo2d

      REAL*8 xmin,xmax,ymin,ymax,zmin,zmax !range dei valori della top.
      REAL*8 ratio

      REAL*8 dist1y,dist2y,dist1x,dist2x,alpha,beta
      INTEGER i,j,k,h,l,ii,jj
      INTEGER ijk
        
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
	   dist1x=cx(i)-xtop(nextx(i)-1)
	   dist2x=xtop(nextx(i))-cx(i)
	   alpha=dist1x/(dist1x+dist2x)
	   DO j=1,ny
	      dist1y=cy(j)-ytop(nexty(j)-1)
	      dist2y=ytop(nexty(j))-cy(j)
	      beta=dist1y/(dist1y+dist2y)

	      topo2d(i,j)=beta*(alpha*ztop2d(nextx(i),nexty(j)) &
		   +(1-alpha)*ztop2d(nextx(i)-1,nexty(j)))      &
		   +(1-beta)*(alpha*ztop2d(nextx(i),nexty(j)-1) &
		   +(1-alpha)*ztop2d(nextx(i)-1,nexty(j)-1))

	   ENDDO
	ENDDO

! cerca i nodi sotto la topografia

	DO i=1,nx
	   DO j=1,ny
	      DO k=1,nz

		IF (cz(k) <= topo2d(i,j)) THEN
		   ord2d(i,j)=k    !ultimo nodo sotto la top.
		ENDIF

                ! ... dist defines implicitly the profile
                !
	        ijk = i + (j-1) * nx + (k-1) * nx * ny
	        dist(ijk) = cz(k) - topo2d(i,j)

	      ENDDO
	   ENDDO
	ENDDO
      
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
