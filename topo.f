!----------------------------------------------------------------------
      MODULE volcano_topography
!
! ... Import the topography from a Digital Elevation Model (DEM)
! ... The topography is discretized into a step-wise geometry and 
! ... the opportune flag is assigned to topographic cells
!
!----------------------------------------------------------------------
      USE dimensions, ONLY: nx, ny, nz, ntot, no
      USE io_files, ONLY: topounit, topofile
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
!
! ... Input parameters
!
      INTEGER :: itp, iavv, ismt
      REAL*8 :: cellsize, filtersize
      REAL*8 :: rim_quota
      REAL *8 :: min_angle, max_angle
      LOGICAL :: nocrater, itrans, seatable, write_improfile
      !
      ! ... file name for the topography
      CHARACTER(LEN=80) :: dem_file
!
! ... Input topography points (from DEM)
!
      REAL*8, ALLOCATABLE :: xtop(:), ytop(:), ztop(:), ztop2d(:,:)
!
! ... Temporary arrays for topographic points
!
      REAL*8, ALLOCATABLE :: xdem(:), ydem(:), zdem(:,:)
!
! ... Coordinates of the vent center
!
      INTEGER :: icenter, jcenter
!
! ... Output: Topographic elevation at the mesh points
! ... (centered and staggered)
!
      REAL*8, ALLOCATABLE  :: topo_c(:)
      REAL*8, ALLOCATABLE  :: topo_x(:)
      REAL*8, ALLOCATABLE  :: topo2d_c(:,:)
      REAL*8, ALLOCATABLE  :: topo2d_x(:,:)
      REAL*8, ALLOCATABLE  :: topo2d_y(:,:)
!
      PUBLIC
      PRIVATE :: icenter,jcenter,xdem,ydem,zdem, dist
      PRIVATE :: dem_header, vdem
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE import_topography
!
! ... Read the topography file
!
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE grid, ONLY: iv, jv

      IMPLICIT NONE
!
      IF (mpime == root) THEN
              OPEN(UNIT=topounit,FILE=topofile,STATUS='UNKNOWN')
              WRITE(topounit,*) 'Report of topographic conditions'
              WRITE(topounit,*)
      END IF
!
! ... Read external topography files.
! ... In 3D translate horizontally the  computational mesh 
! ... (geographic UTM coordinates).
!
      IF (job_type == JOB_TYPE_2D ) THEN

        CALL read_2Dprofile

      ELSE IF (job_type == JOB_TYPE_3D) THEN

        ! ... Read the original topography
        !
        CALL read_dem_ascii

        ! ... Crop the dem and change the resolution + Filtering
        !
        CALL resize_dem
        
        ! ... Radial averaging + Filtering
        !
        IF (iavv >= 1) CALL average_dem

        CALL compute_UTM_coords
        
      END IF
!
! ... Interpolate the topography on the cell centers.
! ... Set the cell flags = 3 below the topographic level.
! ... Compute the distance of each cell center from the ground. 
! ... Translate vertically the mesh.
!
      CALL set_profile
!
      RETURN
      END SUBROUTINE import_topography
!----------------------------------------------------------------------
      SUBROUTINE set_profile
!
      USE array_filters, ONLY: interp
      USE control_flags, ONLY: job_type, lpr
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE grid, ONLY: x, xb, y, yb, z, zb, iv, jv, kv, dzmax
      USE io_files, ONLY: testunit, tempunit
      USE parallel, ONLY: mpime, root

      IMPLICIT NONE
      REAL*8 :: transl_z = 0.D0
      INTEGER :: i,j,k,ijk
      INTEGER :: itopo
!      INTEGER :: ios
      LOGICAL, ALLOCATABLE :: dummy(:)
      REAL*8, ALLOCATABLE  :: topo(:)
      REAL*8, ALLOCATABLE  :: topo2d(:,:)
!
      ALLOCATE (dummy(ntot))
      dummy = .FALSE.
!
      ! ... Allocate and initialize 'dist' (defines implicitly the profile).
      !
      IF (.NOT.ALLOCATED(dist)) ALLOCATE (dist(ntot))
      dist = 1.0D10
!
      IF (job_type == JOB_TYPE_2D ) THEN
        !
        IF (.NOT.ALLOCATED(next)) ALLOCATE(next(nx))
        IF (.NOT.ALLOCATED(ord)) ALLOCATE(ord(nx))
        next = 0; ord = 0
        ALLOCATE(topo(nx))
!
        IF (MAXVAL(xtop) < MAXVAL(x)) &
          CALL error('topo','Computational domain exceeds topography',1)
!
        ! ... Interpolate the topography on the computational mesh 'x'
        ! ... (centered horizontally)
        !
        CALL interp(xtop,ztop,x,topo,next)
!
        DO i = 1, nx
          ! ... Topography is expressed in centimeters
          !
          topo(i) = topo(i) * 1.D2
          itopo = NINT(topo(i))
          topo(i) = itopo * 1.D-2
        ENDDO
!        
        ! ... Translate vertically the mesh to minimize the
        ! ... number of topographic cells
        !
        IF (itrans) transl_z = MINVAL(topo) - zb(1)

        IF( mpime == root ) THEN
          WRITE(topounit,*) 'Translating mesh vertically'
          WRITE(topounit,*) 'Minimum topographic quota: ', MINVAL(topo)
          WRITE(topounit,*) 'Maximum topographic quota: ', MAXVAL(topo)
          WRITE(topounit,*) 'Translation: ', transl_z
          WRITE(topounit,*) 
        END IF
        z  = z  + transl_z
        zb = zb + transl_z
!
        IF (mpime == root) THEN
          OPEN(UNIT=tempunit,FILE='over.log',STATUS='UNKNOWN')
          WRITE(tempunit,*) 'Coordinates of the last cell laying below the topography'
        END IF
        ! ... Re-set the 'ord' array and set the implicit profile
        ! ... 'ord(i)' is the last cell laying below the topography
        !
        DO i = 1, nx
          DO k = 1, nz
            IF (zb(k) <= topo(i)) ord(i) = k  
            ijk = i + (k-1) * nx
            dist(ijk) = z(k) - topo(i)
          END DO
            IF (mpime == root) THEN
                WRITE(tempunit,*) 'i= ', i, 'x(i)=', x(i), 'k= ', ord(i), 'zb(k)= ', zb(ord(i))
            END IF
        END DO
        kv = ord(iv)
!
        DEALLOCATE(topo)
        IF (mpime == root) CLOSE(tempunit)
        !
      ELSE IF (job_type == JOB_TYPE_3D) THEN
        !
        IF (.NOT.ALLOCATED(nextx)) ALLOCATE(nextx(nx))
        IF (.NOT.ALLOCATED(nexty)) ALLOCATE(nexty(ny))
        IF (.NOT.ALLOCATED(ord2d)) ALLOCATE(ord2d(nx,ny))
        nextx = 0; nexty = 0; ord2d = 0
        ALLOCATE(topo2d(nx,ny))
        
        IF (MAXVAL(xtop) < MAXVAL(x)) THEN
          IF( mpime == root ) THEN
            WRITE(topounit,*) MAXVAL(xtop), MAXVAL(x) 
          END IF
          CALL error('topo','Computational domain exceeds topography x',1)
        END IF
        IF (MAXVAL(ytop) < MAXVAL(y)) THEN
          IF( mpime == root ) THEN
            WRITE(topounit,*) MAXVAL(ytop), MAXVAL(y) 
          END IF
          CALL error('topo','Computational domain exceeds topography y',1)
        END IF

        CALL interp(xtop, ytop, ztop2d, x, y, topo2d, nextx, nexty)

        ! ... Topography must be accurate to centimeters
        !
        DO j=1,ny
          DO i=1,nx
            topo2d(i,j) = topo2d(i,j) * 1.D2
            itopo = NINT(topo2d(i,j))
            topo2d(i,j) = itopo * 1.D-2
          END DO
        END DO
        !
        ! ... Translate vertically the numerical mesh to minimize the
        ! ... number of topographic cells
        !
        IF (itrans) transl_z = MINVAL(topo2d) - zb(1)

        IF( mpime == root ) THEN
          WRITE(topounit,*) 'Translating mesh vertically'
          WRITE(topounit,*) 'Minimum topographic quota: ', MINVAL(topo2d)
          WRITE(topounit,*) 'Maximum topographic quota: ', MAXVAL(topo2d)
          WRITE(topounit,*) 'Translation: ', transl_z
          WRITE(topounit,*) 
        END IF
        z  = z  + transl_z
        zb = zb + transl_z
        
        ! ... Reset the 'ord2d' array after translation 
        ! ... and set the implicit profile
        !
        DO j = 1, ny
          DO i = 1, nx
            DO k = 1, nz
              IF (zb(k) <= topo2d(i,j)) ord2d(i,j) = k  
              ijk = i + (j-1) * nx + (k-1) * nx * ny
              dist(ijk) = z(k) - topo2d(i,j)
            END DO
          END DO
        END DO
        kv = ord2d(iv,jv)
        !
        DEALLOCATE(topo2d)
        !
      END IF
!
! ... set cell flags on topography
!
      CALL set_flags

      DEALLOCATE (dummy)

      RETURN
      END SUBROUTINE set_profile
!----------------------------------------------------------------------
      SUBROUTINE read_2Dprofile
      USE parallel, ONLY: mpime, root
      USE io_files, ONLY: tempunit, testunit
      IMPLICIT NONE

      INTEGER :: n, noditop
      CHARACTER(LEN=80) :: topo_file
!
! ... Processor 'root' reads the topography defined on a 
! ... (generic) rectilinear mesh "xtop". Then broadcast
! ... to all processors.
!
      topo_file = TRIM(dem_file)
      IF (mpime == root) THEN
        WRITE(topounit,*) 'Reading DEM file: ', dem_file
        WRITE(topounit,*) 
        OPEN(UNIT=tempunit, FILE=topo_file, STATUS='OLD',ERR=199)
        READ(tempunit,*,ERR=199) noditop
      END IF
!
      CALL bcast_integer(noditop,1,root)
!
      ALLOCATE(xtop(noditop))
      ALLOCATE(ztop(noditop))
!
      IF (mpime == root) THEN
        DO n=1, noditop
          READ(tempunit,*,ERR=199) xtop(n),ztop(n)
          !
          IF (seatable) ztop(n) = MAX(ztop(n),0.D0)
          !
        END DO
        CLOSE(tempunit)
      END IF
!
      CALL bcast_real(xtop,noditop,root)
      CALL bcast_real(ztop,noditop,root)
!
      RETURN
!
 199  CALL error('topo.f', 'error in reading temp unit', tempunit)
!
      END SUBROUTINE read_2Dprofile
!----------------------------------------------------------------------
      SUBROUTINE read_dem_ascii
!
      USE parallel, ONLY: mpime, root
      USE io_files, ONLY: tempunit
      IMPLICIT NONE

      CHARACTER(LEN=80) :: topo_file
      INTEGER :: nodidemx, nodidemy
      ! ... Coordinates of the upper-left (ul) corner
      REAL*8  :: xll, yll, xur, yur, xul, yul, xlr, ylr
      REAL*8  :: dd
      INTEGER :: noval
      INTEGER :: elevation
      INTEGER :: p
      INTEGER :: i, j
      INTEGER :: halfx
!
! ... Processor 'root' reads the standard header of the
! ... ASCII DEM format, then broadcasts it to all processors.
!
      IF (mpime == root) THEN
        topo_file = TRIM(dem_file)
        WRITE(topounit,*) 'Reading topography from file: ', topo_file
        OPEN(UNIT=tempunit, FILE=topo_file, STATUS='OLD',ERR=199)
        READ(tempunit,*,ERR=199) nodidemx
        READ(tempunit,*,ERR=199) nodidemy
        READ(tempunit,*,ERR=199) xll
        READ(tempunit,*,ERR=199) yll
        READ(tempunit,*,ERR=199) dd
        READ(tempunit,*,ERR=199) noval
      END IF
!
      CALL bcast_integer(nodidemx,1,root)
      CALL bcast_integer(nodidemy,1,root)
      CALL bcast_real(xll,1,root)
      CALL bcast_real(yll,1,root)
      CALL bcast_real(dd,1,root)
      CALL bcast_integer(noval,1,root)
!
      vdem%nx           = nodidemx
      vdem%ny           = nodidemy
      vdem%xcorner      = xll
      vdem%ycorner      = yll
      vdem%cellsize     = dd
      vdem%nodata_value = noval
!
! ... If the lower-left coordinates are provided,
! ... use the following...
!
      vdem%xcorner      = xll
      vdem%ycorner      = yll
!
      xul = xll
      yul = yll + vdem%cellsize * (vdem%ny - 1)
      xlr = xll + vdem%cellsize * (vdem%nx - 1)
      xur = xlr
!
! ... These are used to check the domain size
!
!      xur = xul + vdem%cellsize * (vdem%nx - 1)
!      yll = yul - vdem%cellsize * (vdem%ny - 1)
!      xll = xul
!
      ALLOCATE(zdem(vdem%nx,vdem%ny))
      ALLOCATE(xdem(vdem%nx))
      ALLOCATE(ydem(vdem%ny))
!
! ... Compute the UTM coordinates of each element of the topography
! ... from the upper-left corner.
!
      DO i = 1, vdem%nx
        xdem(i) = vdem%xcorner + (i-1) * vdem%cellsize
      END DO
      !
      ! ... upper left corner
      !
      !DO j = vdem%ny, 1, -1
      DO j = 1, vdem%ny
!        ydem(j) = vdem%ycorner - (vdem%ny - j) * vdem%cellsize
        ydem(j) = vdem%ycorner + (j-1) * vdem%cellsize
      END DO
!
      IF (mpime == root) THEN
        WRITE(topounit,*) 'DEM resolution: ', vdem%cellsize, ' [m]'
        WRITE(topounit,*) 'DEM limits: '
        WRITE(topounit,'(2F12.2)') xdem(1), xdem(vdem%nx)
        WRITE(topounit,'(2F12.2)') ydem(1), ydem(vdem%ny)
      END IF
!
! ... Processor 'root' reads the matrix of the elevation
! ... ( in centimeters ) and broadcasts.
!
      IF (mpime == root) THEN
        !DO j = vdem%ny, 1, -1
        !DO j = 1, vdem%ny
          !DO i = 1, vdem%nx
            !READ(tempunit,*,ERR=199) elevation
            READ(tempunit,*,ERR=199) ((zdem(i,j),i=1,vdem%nx),j=vdem%ny, 1, -1)
!            zdem = zdem / 100
            !
            !IF (seatable) zdem(i,j) = MAX(zdem(i,j),0.D0)
            !
          !END DO
        !END DO
        CLOSE(tempunit)
!
!
! ... TEST!!!!
!
!      halfx = vdem%nx/2
!      zdem(i,:)=zdem(i,halfx)
!
      END IF
      CALL bcast_real(zdem,nodidemx*nodidemy,root)
!
      RETURN
!
 199  CALL error('topo.f', 'error in reading temp unit', tempunit)
!
      END SUBROUTINE read_dem_ascii
!----------------------------------------------------------------------
! ... Interpolates or downsizes the DEM with a prescribed uniform 
! ... resolution ('cellsize' is read in input)
!
      SUBROUTINE resize_dem
!
      USE grid, ONLY: domain_x, domain_y, dxmax, dymax, dxmin, dymin
      USE grid, ONLY: center_x, center_y, alpha_x, alpha_y
      USE grid, ONLY: dx, dy
      USE array_filters, ONLY: interp, mean_filter, gaussian_filter, &
                                       barnes_filter
      IMPLICIT NONE
      INTEGER :: noditopx, noditopy
      INTEGER :: nodidemx, nodidemy
      REAL*8 :: newsizex, newsizey
      REAL*8 :: xll, yll, xur, yur, xul, yul
      REAL*8 :: zmin
      REAL*8 :: maxcellsize, mincellsize
      INTEGER, ALLOCATABLE :: ntx(:), nty(:)
      INTEGER :: i,j
!
      ! ... Crop the topography and compute the number of 
      ! ... DEM elements. Allocate elevation arrays.
      !
      newsizex = domain_x + 2.D0 * MAX(dxmax, cellsize)
      newsizey = domain_y + 2.D0 * MAX(dymax, cellsize)
      noditopx = INT(newsizex / cellsize) + 1
      noditopy = INT(newsizey / cellsize) + 1
      !
      ! ... Old dimensions
      nodidemx = vdem%nx
      nodidemy = vdem%ny
!      
      ALLOCATE( xtop(noditopx) )
      ALLOCATE( ytop(noditopy) )
      ALLOCATE( ztop2d(noditopx,noditopy) )
      ALLOCATE( ntx(noditopx) )
      ALLOCATE( nty(noditopy) )
!
      ! ... Compute the UTM coordinates of the lower-left ('x/y-ll'), 
      ! ... the upper-right ('x/y-ur') and the upper-left ('x/y-ul') corners. 
      ! ... Check that the new domain is contained in the input topographic DEM.
      !
      xll = center_x - alpha_x * domain_x - dxmax
      yll = center_y - alpha_y * domain_y - dymax
      !
      xtop(1) = xll
      ytop(1) = yll
      xur = xll + cellsize * (noditopx - 1)
      yur = yll + cellsize * (noditopy - 1)
      xul = xll
      yul = yur
      !
      ! ... Reset the resized DEM parameters as default
      !
      vdem%nx           = noditopx
      vdem%ny           = noditopy
      vdem%xcorner      = xul
      vdem%ycorner      = yul
      vdem%cellsize     = cellsize
!     
      ! ... Compute the new UTM coordinates of each element.
      !
      DO i = 2, vdem%nx
        xtop(i) = xtop(i-1) + vdem%cellsize
      END DO
      DO j = 2, vdem%ny
        ytop(j) = ytop(j-1) + vdem%cellsize
      END DO
!
      IF (mpime == root) THEN
        WRITE(topounit,*) 'DEM is resized'
        WRITE(topounit,*) 'New resolution: ', vdem%cellsize, ' [m]'
        WRITE(topounit,*) 'New DEM limits: '
        WRITE(topounit,'(2F12.2)') xtop(1), xtop(vdem%nx)
        WRITE(topounit,'(2F12.2)') ytop(1), ytop(vdem%ny)
      END IF
!
      IF ( xll < xdem(1) ) &
        CALL error('resize_dem','Domain exceeds West DEM limits',1)
      IF ( xur > xdem(nodidemx) ) &
        CALL error('resize_dem','Domain exceeds East DEM limits',1)
      IF ( yll < ydem(1) ) &
        CALL error('resize_dem','Domain exceeds South DEM limits',1)
      IF ( yul > ydem(nodidemy) ) &
        CALL error('resize_dem','Domain exceeds North DEM limits',1)
!
      ! ... Interpolate the DEM over the new mesh.
      !
      CALL interp(xdem, ydem, zdem, xtop, ytop, ztop2d, ntx, nty)
!
      maxcellsize = MAX(MAXVAL(dx),MAXVAL(dy))
      mincellsize = MIN(MINVAL(dx),MINVAL(dy))
      !filtersize  = MAX(filtersize,0.5D0*maxcellsize)
      !filtersize  = MAX(filtersize,0.5D0*mincellsize)
      IF (filtersize >= cellsize) THEN
         IF (ismt == 1) THEN
                 CALL mean_filter(xtop,ytop,ztop2d,filtersize)
         ELSE IF (ismt == 2) THEN
                 CALL gaussian_filter(xtop,ytop,ztop2d,cellsize,filtersize)
         ELSE IF (ismt == 3) THEN
                 CALL barnes_filter(xtop,ytop,ztop2d,cellsize,filtersize)
         END IF
      END IF

      DEALLOCATE(ntx)
      DEALLOCATE(nty)
      DEALLOCATE(xdem, ydem, zdem)

      RETURN
      END SUBROUTINE resize_dem
!----------------------------------------------------------------------
! ... This routine builds an axisymmetric volcano topography by averaging
! ... the digital elevation model (dem) around the specified x/y-center
!
      SUBROUTINE average_dem
      USE array_filters, ONLY: mean_filter
      USE grid, ONLY: center_x, center_y
      USE io_files, ONLY: tempunit
!
      IMPLICIT NONE
      INTEGER :: icenter, jcenter
      INTEGER :: i, j, distance, m, l
      INTEGER :: dms, counter
      REAL*8, ALLOCATABLE :: av_quota(:)
      REAL*8, ALLOCATABLE :: rad_dist(:), rad_quota(:)
      INTEGER, ALLOCATABLE :: nk(:)
      INTEGER :: noditopx, noditopy
      REAL*8 :: alpha, tanalpha
      REAL*8 :: pi
!
      pi = 4.D0 * DATAN(1.D0)
!
!  Trasformo gli angoli da gradi in radianti
!
      min_angle = min_angle*pi/180.D0
      max_angle = max_angle*pi/180.D0
!
      noditopx = vdem%nx
      noditopy = vdem%ny
!
      DO i = 1, noditopx
        IF (xtop(i) < center_x) icenter = i
      END DO
      DO j = 1, noditopy
        IF (ytop(j) < center_y) jcenter = j
      END DO
!
      dms = noditopx**2 + noditopy**2
      ALLOCATE( av_quota( 0:dms ) )
      ALLOCATE( nk( 0:dms ) )
      av_quota = 0.D0
      nk       = 0
!
      ! ... Average the topography over axisymmetrix layers
      !
      DO j = 1, noditopy
        DO i = 1, noditopx
          IF (j/=jcenter .AND. i/=icenter) THEN
            alpha = ATAN2(REAL(j-jcenter,8),REAL(i-icenter,8))
          ELSE IF (j == jcenter) THEN
            IF (i>icenter) THEN
              alpha = 0.D0
            ELSE IF (i<icenter) THEN
              alpha = pi
            END IF
          ELSE IF (i == icenter) THEN
            IF (j>jcenter) THEN
              alpha = 0.5D0*pi
            ELSE IF (j<jcenter) THEN
              alpha = -0.5D0*pi
            END IF
          END IF
          IF (alpha >= min_angle .AND. alpha <= max_angle) THEN
            distance = (i - icenter)**2 + (j - jcenter)**2
            av_quota(distance) = av_quota(distance) + ztop2d(i,j)
            nk(distance) = nk(distance) + 1
          END IF
        END DO
      END DO
      !
      DO m = 1, dms
        IF( nk(m) > 0 ) av_quota(m) = av_quota(m) / nk(m)
      END DO
      !
      ! ... Count the number of non-zero quotas
      !
      counter = 0
      DO j = 0, dms
        IF (av_quota(j) > 0.D0) counter = counter + 1
      END DO
      IF (counter == 0) RETURN
      !
      ! ... 'rad_dist' is the array of the squared distances of mesh point
      ! ... 'rad_quotas' is the array of averaged quotas at 'rad_dist' locations
      !
      ALLOCATE( rad_dist(counter), rad_quota(counter) )
      rad_dist = 0.D0; rad_quota = 0.D0
      !
      i = 0
      DO l = 0, dms
        IF (av_quota(l) > 0.D0) THEN
                i = i+1
                rad_dist(i) = REAL(l,8)
                rad_quota(i) = av_quota(l)
        END IF
      END DO
!      
      CALL mean_filter(rad_dist, rad_quota, filtersize**2)
      !
      ! ... Write out the 2D averaged profile
      !
      IF (mpime == root) THEN
        OPEN(tempunit,FILE='radial_dem.dat',STATUS='UNKNOWN')
        DO m = 1, SIZE(rad_dist)
            WRITE(tempunit,*) cellsize*DSQRT(rad_dist(m)), rad_quota(m)
        END DO
        CLOSE(tempunit)
      END IF
!
      DO i = 1, counter
        av_quota(NINT(rad_dist(i))) = rad_quota(i)
      END DO
!
      ! ... Assign the new averaged values to the
      ! ... topographic arrays
      !
      DO j = 1, noditopy
        DO i = 1, noditopx
          IF (j/=jcenter .AND. i/=icenter) THEN
            alpha = ATAN2(REAL(j-jcenter,8),REAL(i-icenter,8))
          ELSE IF (j == jcenter) THEN
            IF (i>icenter) THEN
              alpha = 0.D0
            ELSE IF (i<icenter) THEN
              alpha = pi
            END IF
          ELSE IF (i == icenter) THEN
            IF (j>jcenter) THEN
              alpha = 0.5D0*pi
            ELSE IF (j<jcenter) THEN
              alpha = -0.5D0*pi
            END IF
          END IF
          IF (alpha >= min_angle .AND. alpha <= max_angle) THEN
            distance = (i - icenter)**2 + (j - jcenter)**2
            ztop2d(i,j) = av_quota(distance)
          END IF
        END DO
      END DO
!
      DEALLOCATE(av_quota)
      DEALLOCATE(rad_quota, rad_dist)
      DEALLOCATE(nk)
!
      RETURN
      END SUBROUTINE average_dem
!----------------------------------------------------------------------
! ... Translates the computational mesh accordingly to the
! ... specified UTM coordinates of the DEM
!
      SUBROUTINE compute_UTM_coords

      USE grid, ONLY: center_x, center_y
      USE grid, ONLY: x, y, xb, yb
      USE grid, ONLY: iv, jv
      IMPLICIT NONE
      
      REAL*8 :: transl_x, transl_y
!
      transl_x = center_x - x(iv)
      transl_y = center_y - y(jv)
!
! ... Translate all mesh arrays horizontally
!
      x  = x  + transl_x
      xb = xb + transl_x
      y  = y  + transl_y
      yb = yb + transl_y
!
      RETURN
      END SUBROUTINE compute_UTM_coords
!----------------------------------------------------------------------
      SUBROUTINE flatten_dem_vent(xvent,yvent,base_radius,crater_radius,quota)
      USE grid, ONLY: zb
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: xvent,yvent,base_radius,crater_radius
      INTEGER, INTENT(IN) :: quota
      REAL*8 :: distance2
      INTEGER, ALLOCATABLE :: rim_east(:), rim_west(:)
      INTEGER, ALLOCATABLE :: rim_north(:), rim_south(:)
      INTEGER :: i,j,k
!
      IF (nocrater) THEN
        !
        ALLOCATE( rim_east(vdem%ny), rim_west(vdem%ny) )
        ALLOCATE( rim_north(vdem%nx), rim_south(vdem%nx) )
        !
        DO j = 1, vdem%ny
          rim_west(j) = vdem%nx 
          rim_east(j) = 1
          DO i = 1, vdem%nx
            distance2 = (xtop(i)-xvent)**2 + (ytop(j)-yvent)**2 
            IF ( (distance2 < crater_radius**2)   .AND. &
                 (ztop2d(i,j) > zb(quota) ) ) THEN
                rim_west(j) = MIN(i,rim_west(j))
                rim_east(j) = MAX(i,rim_east(j))
            END IF
          END DO
        END DO
        !
        DO i = 1, vdem%nx
          rim_south(i) = vdem%ny 
          rim_north(i) = 1
          DO j = 1, vdem%ny
            distance2 = (xtop(i)-xvent)**2 + (ytop(j)-yvent)**2 
            IF ( (distance2 < crater_radius**2)   .AND. &
                 (ztop2d(i,j) > zb(quota) ) ) THEN
                rim_south(i) = MIN(j,rim_south(i))
                rim_north(i) = MAX(j,rim_north(i))
            END IF
          END DO
        END DO
        !
        DO j = 1, vdem%ny
          DO i = 1, vdem%nx
            IF( i >= rim_west(j) .AND. i <= rim_east(j) .AND. &
                j >= rim_south(i) .AND. j <= rim_north(i) ) THEN
                ztop2d(i,j) = zb(quota)
            END IF
          END DO
        END DO
        !
        DEALLOCATE( rim_east, rim_west )
        DEALLOCATE( rim_north, rim_south )
        !
      ELSE
        !
        DO j = 1, vdem%ny
          DO i = 1, vdem%nx
            distance2 = (xtop(i)-xvent)**2 + (ytop(j)-yvent)**2 
            IF( distance2 < base_radius**2 ) THEN
              ztop2d(i,j) = zb(quota)
            END IF
          END DO
        END DO
!
      END IF
!
      ! ... Interpolate the topography on the cell centers.
      ! ... Re-set the cell flags at the base of the crater
      ! ... and the 'ord2d' and 'dist' arrays
      !
      CALL set_profile
!
      RETURN
      END SUBROUTINE flatten_dem_vent
!----------------------------------------------------------------------
      SUBROUTINE flatten_dem_vent_2D(xvent,vent_radius,quota)
      USE grid, ONLY: zb
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: xvent,vent_radius
      INTEGER, INTENT(IN) :: quota
      REAL*8 :: distance2
      INTEGER :: i,k
!
! ... All cells included in the vent radius (plus one)
! ... have the same topographic quota
!
      DO i = 1, SIZE(xtop)
        distance2 = (xtop(i)-xvent)**2
        IF( distance2 < vent_radius**2 ) THEN
          ztop(i) = zb(quota)
          ztop(i+1) = zb(quota)
        END IF
      END DO
!
      ! ... Interpolate the topography on the cell centers.
      ! ... Re-set the cell flags at the base of the crater
      ! ... and the 'ord2d' and 'dist' arrays
      !
      CALL set_profile
!
      RETURN
      END SUBROUTINE flatten_dem_vent_2D
!----------------------------------------------------------------------
      SUBROUTINE flatten_dem_dome(xdome,ydome,dome_radius,quota)
      USE grid, ONLY: zb
      USE io_files, ONLY: testunit, tempunit
      USE parallel, ONLY: mpime, root
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: xdome,ydome,dome_radius
      INTEGER, INTENT(IN) :: quota
      REAL*8 :: distance2, ztp
      REAL*8 :: modified_xdome, modified_ydome
      INTEGER :: i,j,k
!
      modified_xdome = xdome
      modified_ydome = ydome + 0.5*dome_radius
      DO j = 1, vdem%ny
        DO i = 1, vdem%nx
          distance2 = (xtop(i)-modified_xdome)**2 + (ytop(j)-modified_ydome)**2 
          IF( distance2 < dome_radius**2) THEN
            ztp = ztop2d(i,j)
            ztop2d(i,j) = MIN(zb(quota),ztp)
          END IF
        END DO
      END DO
      IF (mpime==root) THEN
        OPEN(tempunit,FILE='DEM_modified.dat',STATUS='UNKNOWN')
       
        WRITE(tempunit,*) vdem%nx    
        WRITE(tempunit,*) vdem%ny    
        WRITE(tempunit,*) vdem%xcorner
        WRITE(tempunit,*) vdem%ycorner
        WRITE(tempunit,*) vdem%cellsize
        WRITE(tempunit,*) '-9999'
        DO j = vdem%ny, 1, -1
          DO i = 1, vdem%nx
            WRITE(tempunit,*) ztop2d(i,j)
          END DO
        END DO
        CLOSE(tempunit)
      END IF
!
      ! ... Interpolate the topography on the cell centers.
      ! ... Re-set the cell flags at the base of the crater
      ! ... and the 'ord2d' and 'dist' arrays
      !
      CALL set_profile
!
      RETURN
      END SUBROUTINE flatten_dem_dome
!----------------------------------------------------------------------
      SUBROUTINE interpolate_profile(cx, cz, topo, ff)
      USE array_filters, ONLY: interp
!
! ... interpolate the topographic profile on a given
! ... computational mesh (either centered or staggered)
! ... as defined by the arrays (cx, cy, cz)

      IMPLICIT NONE

      REAL*8, INTENT(IN), DIMENSION(:) :: cx, cz
      LOGICAL, INTENT(OUT), DIMENSION(:) :: ff
      REAL*8, INTENT(OUT), DIMENSION(:) :: topo

      INTEGER :: i, k, ijk, n
      INTEGER :: itopo
      REAL*8 :: grad
!
! ... locate the grid points near the topographic points
! ... and interpolate linearly the profile  
!
      ff = .FALSE.
!
! ... interpolate the topography on the mesh
!
      CALL interp(xtop,ztop,cx,topo,next)

      DO i = 1, nx
        ! ... Topography is expressed in centimeters
        !
        topo(i) = topo(i) * 1.D2
        itopo = NINT(topo(i))
        topo(i) = itopo * 1.D-2
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
      DO i = 2, nx - 1
        IF ((ord(i)-ord(i-1) > 0).AND.(topo(i+1)-topo(i-1) >= 0)) THEN
          DO k = ord(i-1), ord(i)
            ijk = i + (k-1) * nx
            ff(ijk) = .TRUE.
          END DO
        ELSEIF (ord(i+1)-ord(i) < 0) THEN
          DO k = ord(i+1), ord(i)
            ijk = i + (k-1) * nx
            ff(ijk) = .TRUE.
          END DO
        ELSE
          k = ord(i)
          ijk = i + (k-1) * nx
          ff(ijk) = .TRUE.
        END IF
      END DO
!
      RETURN
      END SUBROUTINE interpolate_profile
!----------------------------------------------------------------------
      SUBROUTINE interpolate_dem(cx, cy, cz, topo2d, ff)
      USE dimensions, ONLY: nx, ny, nz
      USE array_filters, ONLY: interp
      IMPLICIT NONE

      REAL*8, INTENT(IN), DIMENSION(:) :: cx, cy, cz
      REAL*8, INTENT(OUT), DIMENSION(:,:) :: topo2d
      LOGICAL, INTENT(OUT), DIMENSION(:) :: ff

      INTEGER :: i,j,k,ijk
      INTEGER :: itopo

      CALL interp(xtop, ytop, ztop2d, cx, cy, topo2d, nextx, nexty)

      ! ... Topography must be accurate to centimeters
      !
!      DO j=1,ny
!        DO i=1,nx
!          topo2d(i,j) = topo2d(i,j) * 1.D2
!          itopo = NINT(topo2d(i,j))
!          topo2d(i,j) = itopo * 1.D-2
!        END DO
!      END DO

! ... Locate cells laying below the topography 'ord2d(i,j)'

      DO j=1,ny
         DO i=1,nx
            DO k=1,nz
              IF (cz(k) <= topo2d(i,j)) ord2d(i,j) = k  
            ENDDO
         ENDDO
      ENDDO
!
! ... Identify forcing points
! ... (skip boundaries)
!
      DO i = 2, nx - 1
        DO j = 2, ny - 1
          DO k = 2, ord2d(i,j) - 1
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
      SUBROUTINE set_flags
!
! ... Set cell-flag = 3 in cells laying below the topography
! ... "ord" and "ord2d" are the last cell COMPLETELY below the topography
!
      USE control_flags, ONLY: lpr, job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE grid, ONLY: fl, zb, z
      USE grid, ONLY: inlet_cell, noslip_wall, fluid, bl_cell
      USE io_files, ONLY: testunit
      USE flux_limiters, ONLY: muscl
      IMPLICIT NONE

      INTEGER :: i, j, k, ijk, ii, jj
      INTEGER :: q
!
! ... 2D
! ... Define the topography through the no-slip cells
!
      IF( job_type == JOB_TYPE_2D ) THEN
        DO i = 1, nx
          q = ord(i)
          DO k = 1, nz
            ijk = i + (k-1) * nx
            IF (q >= k) THEN
                    IF( fl(ijk) /= inlet_cell ) fl(ijk) = noslip_wall
            ELSE
                    IF( fl(ijk) == noslip_wall ) fl(ijk) = fluid
            END IF
          END DO
        END DO
!
! ... Define a 'boundary layer' of thickness == 3 above the topography
!
        DO i = 1, nx
          q = ord(i)
          DO k = q+1, q+3
            DO ii=-3,+3
              ijk = (i+ii) + (k-1) * nx
                  IF( fl(ijk) == fluid ) fl(ijk) = bl_cell
            END DO
          END DO
        END DO
!
! ... 3D
! ... Define the topography through the no-slip cells
!
      ELSE IF( job_type == JOB_TYPE_3D) THEN
        DO j = 1, ny
          DO i = 1, nx
            q = ord2d(i,j)
            DO k = 1, nz
              ijk = i + (j-1) * nx + (k-1) * nx * ny
              IF (k <= q) THEN
                      IF( fl(ijk) /= inlet_cell ) fl(ijk) = noslip_wall
              ELSE
                      IF( fl(ijk) == noslip_wall ) fl(ijk) = fluid
              END IF
            END DO
          END DO
        END DO
!
! ... Define a layer of thickness == 3 above the topography
!
        DO j = 1, ny
          DO i = 1, nx
            q = ord2d(i,j)
            DO k = q+1, q+3
              DO jj=-3,+3
                DO ii=-3,+3
                  ijk = (i+ii) + (j+jj-1) * nx + (k-1) * nx * ny
                      IF( fl(ijk) == fluid ) fl(ijk) = bl_cell
                END DO
              END DO
            END DO
          END DO
        END DO
!
      END IF
!
      RETURN
      END SUBROUTINE set_flags
!----------------------------------------------------------------------
      SUBROUTINE export_topography
!
! ... Store the value of the discrete topography
! ... If Immersed Boundaries are used, these values
! ... are modified later with more accurate interpolations
!
      USE control_flags, ONLY: job_type, lpr
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE grid, ONLY: fl, noslip_wall, zb, x, bl_cell, fluid, inlet_cell
      USE io_files, ONLY: tempunit, testunit
      IMPLICIT NONE

      INTEGER :: i, j, k, ijk, q
!
      IF( job_type == JOB_TYPE_2D ) THEN
        ALLOCATE(topo_c(nx))
        ALLOCATE(topo_x(nx))
        DO i = 1, nx
          DO k = 1, nz
            ijk = i + (k-1) * nx
            IF( fl(ijk) == noslip_wall) THEN
              topo_c(i) = zb(k)
              q = k
            ELSE IF (fl(ijk) == fluid) THEN
              IF (k <= q + 3) fl(ijk) = bl_cell
            END IF
          END DO
          topo_x(i) = topo_c(i)
        END DO
      ELSE IF( job_type == JOB_TYPE_3D) THEN
        ALLOCATE(topo2d_c(nx,ny))
        ALLOCATE(topo2d_x(nx,ny))
        ALLOCATE(topo2d_y(nx,ny))
        DO j = 1, ny
          DO i = 1, nx
            DO k = 1, nz
              ijk = i + (j-1) * nx + (k-1) * nx * ny
              IF( fl(ijk) == noslip_wall) THEN
                topo2d_c(i,j) = zb(k)
                q = k
              ELSE IF (fl(ijk) == fluid) THEN
                IF (k <= q + 3) fl(ijk) = bl_cell
              END IF
            END DO
            topo2d_x(i,j) = topo2d_c(i,j)
            topo2d_y(i,j) = topo2d_c(i,j)
          END DO
        END DO
      END IF
      !
      ! ... Write a new DEM file (z(x) in 2D, z(x,y) in 3D) for
      ! ... visualization of the topography
      !
      IF (mpime == root .AND. itp>0) THEN
        OPEN(tempunit, FILE='export_topography.dat',STATUS='UNKNOWN')
        WRITE(topounit,*) 'The new DEM is written in file "export_topography.dat"'
        IF( job_type == JOB_TYPE_2D ) THEN
                  DO i=1,nx
                    WRITE(tempunit,'(2F15.6)') x(i), topo_c(i)
                  END DO
        ELSE IF( job_type == JOB_TYPE_3D) THEN
                  WRITE(tempunit,'(10F15.6)') topo2d_c(:,:)
        END IF
        CLOSE(tempunit)
      END IF
!
      RETURN
      END SUBROUTINE export_topography
!----------------------------------------------------------------------
      SUBROUTINE write_profile

      USE control_flags, ONLY: job_type, lpr
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions
      USE grid, ONLY: x, xb, y, yb, z, zb
      USE grid, ONLY: iob, fl, noslip_wall, filled_cell_1, filled_cell_2
      USE io_files, ONLY: tempunit
!
      IMPLICIT NONE
      INTEGER :: n, i, j, k, ijk, cntz
!
! ... Write out the georeferenced mesh coordinates
!
      IF (mpime == root) THEN
      OPEN(tempunit,FILE='mesh.dat')
      WRITE(tempunit,*) 'Georeferenced x-y mesh'
      WRITE(tempunit,*) 'x'
      WRITE(tempunit,14) x
      WRITE(tempunit,*) 'xb'
      WRITE(tempunit,14) xb
      WRITE(tempunit,*) 'y'
      WRITE(tempunit,14) y
      WRITE(tempunit,*) 'yb'
      WRITE(tempunit,14) yb
      WRITE(tempunit,*) 'z'
      WRITE(tempunit,14) z
      WRITE(tempunit,*) 'zb'
      WRITE(tempunit,14) zb
      CLOSE(tempunit)
      END IF
!
 14   FORMAT(5(F20.6))
!
! ... Write out the implicit profile
!
      IF (mpime == root .AND. write_improfile) THEN
        OPEN(UNIT=tempunit,FILE='improfile.dat',FORM='UNFORMATTED')
!        WRITE(tempunit,fmt='(F9.2)') dist
        WRITE(tempunit) dist
        CLOSE(tempunit)
      END IF
!
! ... Control that the cells below the specified_flow blocks
! ... are noslip cells (they could have been modified by the 
! ... immersed boundary procedure)
!
      IF( job_type == JOB_TYPE_2D ) THEN
        DO n = 1, no
          IF (iob(n)%typ == 5) THEN
            DO i = iob(n)%xlo, iob(n)%xhi
              DO k = 1, iob(n)%zlo - 1
                ijk = i + (k-1) * nx
                fl(ijk) = noslip_wall
              END DO
            END DO
          END IF
        END DO
      ELSE IF (job_type == JOB_TYPE_3D) THEN
        DO n = 1, no
          IF (iob(n)%typ == 5) THEN
            DO i = iob(n)%xlo, iob(n)%xhi
              DO j = iob(n)%ylo, iob(n)%yhi
                DO k = 1, iob(n)%zlo - 1
                  ijk = i + (j-1) * nx + (k-1) * nx * ny
                  fl(ijk) = noslip_wall
                END DO
              END DO
            END DO
          END IF
        END DO
      END IF
!
! ... Deallocate all arrays in the topography module
!
      IF (job_type == JOB_TYPE_2D ) THEN
        DEALLOCATE (next)
        DEALLOCATE (ord)
        DEALLOCATE (dist)
        DEALLOCATE (xtop, ztop)
      ELSE IF (job_type == JOB_TYPE_3D) THEN
        DEALLOCATE(nextx, nexty)
        DEALLOCATE(ord2d)
        DEALLOCATE (dist)
        DEALLOCATE (xtop, ytop, ztop2d)
      END IF
!
      RETURN
      END SUBROUTINE write_profile
!----------------------------------------------------------------------
      END MODULE volcano_topography
!----------------------------------------------------------------------
