!----------------------------------------------------------------------
      MODULE volcano_topography
!----------------------------------------------------------------------
      USE dimensions, ONLY: nx, ny, nz, ntot, no
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

! ... input topography points
!
      REAL*8, ALLOCATABLE :: xtop(:), ytop(:), ztop(:), ztop2d(:,:)
      REAL*8, ALLOCATABLE :: xdem(:), ydem(:), zdem(:,:)
      INTEGER :: icenter, jcenter
!
      ! ... Input parameters
      INTEGER :: itp, iavv
      REAL*8 :: cellsize, filtersize
      REAL*8 :: rim_quota
      LOGICAL :: flatten_crater
!
! ... file name for the topography
      CHARACTER(LEN=80) :: dem_file
!
      PUBLIC
      PRIVATE :: icenter,jcenter,xdem,ydem,zdem
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE import_topography
!
! ... Read the topography file
!
      USE control_flags, ONLY: job_type
      USE grid, ONLY: fl
      USE parallel, ONLY: mpime, root

      IMPLICIT NONE

      IF( mpime == root ) THEN
        WRITE(6,*)
        WRITE(6,*) 'Importing topographic matrix ...'
      END IF
!
! ... Read external topography files.
! ... In 3D translate horizontally the  computational mesh 
! ... (geographic UTM coordinates).
!
      IF (job_type == '2D') THEN

        CALL read_2Dprofile

      ELSE IF (job_type == '3D') THEN
        ! ... Read the original topography
        !
        CALL read_dem_ascii

        ! ... Crop the dem and change the resolution
        !
        CALL resize_dem
        
        IF (iavv >= 1) THEN
                ! ... Radial averaging + Filtering
                !
                CALL average_dem
        ELSE
                ! ... Just Filter the high frequencies
                !
                CALL filter_2d(xtop,ytop,ztop2d,100,100)
        END IF

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
      IF( mpime == root ) THEN
        WRITE(6,*) 'END Import topography'
      END IF
!
      RETURN
      END SUBROUTINE import_topography
!----------------------------------------------------------------------
      SUBROUTINE read_2Dprofile
      USE parallel, ONLY: mpime, root
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
        OPEN(UNIT=3, FILE=topo_file, STATUS='OLD')
        READ(3,*) noditop
      END IF
!
      CALL bcast_integer(noditop,1,root)
!
      ALLOCATE(xtop(noditop))
      ALLOCATE(ztop(noditop))
!
      IF (mpime == root) THEN
        DO n=1, noditop
          READ(3,*) xtop(n),ztop(n)
        END DO
        CLOSE(3)
      END IF
!
      CALL bcast_real(xtop,noditop,root)
      CALL bcast_real(ztop,noditop,root)
!
      RETURN
      END SUBROUTINE read_2Dprofile
!----------------------------------------------------------------------
      SUBROUTINE read_dem_ascii
!
      USE parallel, ONLY: mpime, root
      IMPLICIT NONE

      CHARACTER(LEN=80) :: topo_file
      INTEGER :: nodidemx, nodidemy
      ! ... Coordinates of the upper-left (ul) corner
      REAL*8  :: xul, yul
      REAL*8  :: dd
      INTEGER :: noval
      INTEGER :: elevation
      INTEGER :: p
      INTEGER :: i, j
!
! ... Processor 'root' reads the standard header of the
! ... ASCII DEM format, then broadcasts it to all processors.
!
      IF (mpime == root) THEN
        topo_file = TRIM(dem_file)
        OPEN(UNIT=3, FILE=topo_file, STATUS='OLD')
        WRITE(6,*) 'Reading topography file: ', topo_file
        READ(3,*) nodidemx
        READ(3,*) nodidemy
        READ(3,*) xul
        READ(3,*) yul
        READ(3,*) dd
        READ(3,*) noval
      END IF
!
      CALL bcast_integer(nodidemx,1,root)
      CALL bcast_integer(nodidemy,1,root)
      CALL bcast_real(xul,1,root)
      CALL bcast_real(yul,1,root)
      CALL bcast_real(dd,1,root)
      CALL bcast_integer(noval,1,root)
!
      vdem%nx           = nodidemx
      vdem%ny           = nodidemy
      vdem%xcorner      = xul
      vdem%ycorner      = yul
      vdem%cellsize     = dd
      vdem%nodata_value = noval
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
      DO j = vdem%ny, 1, -1
        ydem(j) = vdem%ycorner - (vdem%ny - j) * vdem%cellsize
      END DO
!
! ... Processor 'root' reads the matrix of the elevation
! ... ( in centimeters ) and broadcasts.
!
      IF (mpime == root) THEN
        DO j = vdem%ny, 1, -1
          DO i = 1, vdem%nx
            READ(3,*) elevation
            zdem(i,j) = DBLE(elevation) / 100.D0
          END DO
        END DO
      END IF
!
      CALL bcast_real(zdem,nodidemx*nodidemy,root)
!
      RETURN
      END SUBROUTINE read_dem_ascii
!----------------------------------------------------------------------
! ... Interpolates or downsizes the DEM with a prescribed uniform 
! ... resolution ('cellsize' is read in input)
!
      SUBROUTINE resize_dem
!
      USE grid, ONLY: domain_x, domain_y, dxmax, dymax, dxmin, dymin
      USE grid, ONLY: center_x, center_y, alpha_x, alpha_y
      IMPLICIT NONE
      INTEGER :: noditopx, noditopy
      INTEGER :: nodidemx, nodidemy
      REAL*8 :: newsizex, newsizey
      REAL*8 :: xll, yll, xur, yur, xul, yul
      INTEGER, ALLOCATABLE :: ntx(:), nty(:)
      INTEGER :: i,j
!
      ! ... Crop the topography and compute the number of 
      ! ... DEM elements. Allocate elevation arrays.
      !
      newsizex = domain_x + 2.D0 * dxmax
      newsizey = domain_y + 2.D0 * dymax
      noditopx = INT(newsizex / cellsize) + 1
      noditopy = INT(newsizey / cellsize) + 1
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
      xtop(1) = MAX(xdem(1),xll)
      ytop(1) = MAX(ydem(1),yll)
      xur = xll + cellsize * (noditopx - 1)
      yur = yll + cellsize * (noditopy - 1)
      xul = xll
      yul = yur
      nodidemx = vdem%nx
      nodidemy = vdem%ny
      IF (xur > xdem(nodidemx)) &
        noditopx = INT( (xdem(nodidemx) - xtop(1)) / cellsize) 
      IF (yur > ydem(nodidemy)) &
        noditopy = INT( (ydem(nodidemy) - ytop(1)) / cellsize) 
!     
      ! ... Reset the resized DEM parameters as default
      !
      IF (mpime == root) THEN
        WRITE(6,*) 'DEM is resized'
        WRITE(6,*) 'Old resolution: ', vdem%cellsize, ' [m]'
      END IF
!
      vdem%nx           = noditopx
      vdem%ny           = noditopy
      vdem%xcorner      = xul
      vdem%ycorner      = yul
      vdem%cellsize     = cellsize
!
      IF (mpime == root) THEN
        WRITE(6,*) 'New resolution: ', vdem%cellsize, ' [m]'
      END IF
!
      ! ... Compute the new UTM coordinates of each element.
      !
      DO i = 2, noditopx
        xtop(i) = xtop(i-1) + cellsize
      END DO
      DO j = 2, noditopy
        ytop(j) = ytop(j-1) + cellsize
      END DO
!
      ! ... Interpolate the DEM over the new mesh.
      !
      CALL interp_2d(xdem, ydem, zdem, xtop, ytop, ztop2d, ntx, nty)
!
      ! ... Find the closest topographic element
      ! ... to the prescribed vent center.
      !
      DO i = 1, noditopx
        IF (xtop(i) <= center_x) icenter = i
      END DO
      DO j = 1, noditopy
        IF (ytop(j) <= center_y) jcenter = j
      END DO
!      
      ! ... The vent center must coincide with the
      ! ... topographic element.
      !
      center_x = xtop(icenter)
      center_y = ytop(jcenter)
!
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
!
      IMPLICIT NONE
      INTEGER :: i, j, distance, m, l
      INTEGER :: dms, samp, counter
      REAL*8, ALLOCATABLE :: av_quota(:)
      REAL*8, ALLOCATABLE :: rad_dist(:), rad_quota(:)
      INTEGER, ALLOCATABLE :: nk(:)
      INTEGER :: noditopx, noditopy
!
      noditopx = vdem%nx
      noditopy = vdem%ny
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
          distance = (i - icenter)**2 + (j - jcenter)**2
          av_quota(distance) = av_quota(distance) + ztop2d(i,j)
          nk(distance) = nk(distance) + 1
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
      samp = 100
      !
      ! ... 'rad_dist' is the sorted array of the squared distances of mesh point
      ! ... 'rad_quotas' is the array of averaged quotas at 'rad_dist' locations
      !
      ALLOCATE( rad_dist(counter), rad_quota(counter) )
      rad_dist = 0.D0; rad_quota = 0.D0
      !
      i = 0
      DO j = 0, dms
        IF (av_quota(j) > 0.D0) THEN
                i = i+1
                rad_dist(i) = REAL(j,8)
                rad_quota(i) = av_quota(j)
        END IF
      END DO
!      
      CALL filter_1d(rad_dist, rad_quota, samp)

      DO i = 1, counter
        av_quota(NINT(rad_dist(i))) = rad_quota(i)
      END DO
!
      ! ... Assign the new averaged values to the
      ! ... topographic arrays
      !
      DO j = 1, noditopy
        DO i = 1, noditopx
          distance = (i - icenter)**2 + (j - jcenter)**2
          ztop2d(i,j) = av_quota(distance)
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
      SUBROUTINE flatten_dem(xvent,yvent,base_radius,crater_radius,quota)
      USE grid, ONLY: zb
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: xvent,yvent,base_radius,crater_radius
      INTEGER, INTENT(IN) :: quota
      REAL*8 :: distance2
      INTEGER, ALLOCATABLE :: rim_east(:), rim_west(:)
      INTEGER, ALLOCATABLE :: rim_north(:), rim_south(:)
      INTEGER :: i,j,k
!
      IF (flatten_crater) THEN
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
      ! ... Re-set the cell flags at the base of the crater
      ! ... and the 'ord2d' and 'dist' arrays
      !
      DEALLOCATE(nextx,nexty,dist,ord2d)
      CALL set_profile
!
! ... Write out the new DEM file
      OPEN(17,FILE='newdem.dat')
      WRITE(17,*) vdem%nx
      WRITE(17,*) vdem%ny
      WRITE(17,*) vdem%xcorner
      WRITE(17,*) vdem%ycorner
      WRITE(17,*) vdem%cellsize
      WRITE(17,*) vdem%nodata_value
      WRITE(17,*) NINT(ztop2d*100.D0)
      CLOSE(17)
!
      RETURN
      END SUBROUTINE flatten_dem
!----------------------------------------------------------------------
! ... Filter out high frequency modes by successively subsampling,
! ... averaging and interpolating 
!
      SUBROUTINE filter_1d(x1,f1,nsm)

      IMPLICIT NONE
      REAL*8, INTENT(INOUT), DIMENSION(:) :: x1, f1
      INTEGER, INTENT(IN) :: nsm
      REAL*8, ALLOCATABLE, DIMENSION(:) :: x2, f2
      INTEGER, ALLOCATABLE :: dummy(:)
      REAL*8 :: sstep
      INTEGER :: i,j
      INTEGER :: counter, cnt

      counter = SIZE(x1)

      ! ... 'x2' is a uniform subsample of 'x1'
      ! ... 'f2' is an average of 'f1' at 'x2' locations 
      ALLOCATE( x2(nsm), f2(nsm) )
      ALLOCATE( dummy(counter) )

      x2 = 0.D0; f2 = 0.D0
      ! ... Uniformly Subsample and smooth the radial function 
      ! ... of the averaged quota
      !
      sstep = REAL(x1(counter)/(nsm-1),8)
      x2(1) = x1(1)
      f2(1) = f1(1)
      DO i = 2, nsm-1
        x2(i) = x2(i-1) + sstep
        cnt = 0
        jloop: DO j = 1, counter
          IF (DABS(x1(j)-x2(i)) <= sstep ) THEN
                  f2(i) = f2(i) + f1(j)
                  cnt = cnt + 1
          END IF
        END DO jloop
        f2(i) = f2(i) / cnt
      END DO
      x2(nsm) = x1(counter)
      f2(nsm) = f1(counter)
!
      ! ... Linearly interpolate quotas
      !
      CALL interp_1d(x2, f2, x1, f1, dummy)

      DEALLOCATE(x2, f2, dummy)
      RETURN
      END SUBROUTINE filter_1d
!----------------------------------------------------------------------
! ... Filter out high frequency modes by successively subsampling,
! ... averaging and interpolating 
!
      SUBROUTINE filter_2d(x1,y1,f1,nsmx,nsmy)

      IMPLICIT NONE
      REAL*8, INTENT(INOUT), DIMENSION(:) :: x1, y1
      REAL*8, INTENT(INOUT), DIMENSION(:,:) :: f1
      INTEGER, INTENT(IN) :: nsmx, nsmy
      REAL*8, ALLOCATABLE, DIMENSION(:) :: x2, y2
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: f2
      INTEGER, ALLOCATABLE :: dmx(:), dmy(:)
      INTEGER, ALLOCATABLE :: dcx(:,:), dcy(:,:)
      REAL*8 :: sstepx, sstepy
      INTEGER :: i,j,cx,cy,c
      INTEGER :: counterx, countery, cnt

      counterx = SIZE(x1)
      countery = SIZE(y1)

      ! ... 'x2' is a uniform subsample of 'x1'
      ! ... 'f2' is an average of 'f1' at 'x2' locations 
      ALLOCATE( x2(nsmx), y2(nsmy), f2(nsmx,nsmy) )
      ALLOCATE( dmx(counterx), dmy(countery)  )
      ALLOCATE( dcx(counterx,2), dcy(countery,2)  )

      x2 = 0.D0; y2 = 0.D0; f2 = 0.D0
      dcx(:,1) = counterx; dcy(:,1) = countery 
      dcx(:,2) = 0; dcy(:,2) = 0
      ! ... Uniformly Subsample and smooth the x-y function 
      ! ... of the quota
      !
      sstepx = REAL( (x1(counterx)-x1(1))/(nsmx-1) , 8 )
      sstepy = REAL( (y1(countery)-y1(1))/(nsmy-1) , 8 )

      ! ... corners
      !
      ! ... Coordinates of the subsampled set
      !
      x2(1) = x1(1)
      y2(1) = y1(1)
      !
      DO i = 2, nsmx-1
        x2(i) = x2(i-1) + sstepx
      END DO
      DO j = 2, nsmy-1
        y2(j) = y2(j-1) + sstepy
      END DO
      !
      x2(nsmx) = x1(counterx)
      y2(nsmy) = y1(countery)
!      
      f2(1,1) = f1(1,1)
      f2(nsmx,nsmy) = f1(counterx,countery)
      f2(1,nsmy) = f1(1,countery)
      f2(nsmx,1) = f1(counterx,1)

      ! ... 1D average on boundaries
      !
      DO i = 2, nsmx-1
        cnt = 0
        DO c = 1, counterx
          IF (DABS(x1(c)-x2(i)) <= sstepx ) THEN
                  f2(i,1) = f2(i,1) + f1(c,1)
                  f2(i,nsmy) = f2(i,nsmy) + f1(c,countery)
                  cnt = cnt + 1
                  dcx(i,1) = MIN(c,dcx(i,1))
                  dcx(i,2) = MAX(c,dcx(i,2))
          END IF
        END DO
        f2(i,1) = f2(i,1) / cnt
        f2(i,nsmy) = f2(i,nsmy) / cnt
      END DO
      !
      DO j = 2, nsmy-1
        cnt = 0
        DO c = 1, countery
          IF (DABS(y1(c)-y2(j)) <= sstepy ) THEN
                  f2(1,j) = f2(1,j) + f1(1,c)
                  f2(nsmx,j) = f2(nsmx,j) + f1(counterx,c)
                  cnt = cnt + 1
                  dcy(j,1) = MIN(c,dcy(j,1))
                  dcy(j,2) = MAX(c,dcy(j,2))
          END IF
        END DO
        f2(1,j) = f2(1,j) / cnt
        f2(nsmx,j) = f2(nsmx,j) / cnt
      END DO
      !
      DO j = 2, nsmy-1
        DO i = 2, nsmx-1
          cnt = 0
          DO cy = dcy(j,1), dcy(j,2)
            DO cx = dcx(i,1), dcx(i,2)
              IF (DABS(y1(cy)-y2(j)) <= sstepy   .AND. &
                  DABS(x1(cx)-x2(i)) <= sstepx ) THEN
                      f2(i,j) = f2(i,j) + f1(cx,cy)
                      cnt = cnt + 1
              END IF
            END DO
          END DO
          f2(i,j) = f2(i,j) / cnt
        END DO
      END DO
!
      ! ... Linearly interpolate quotas
      !
      CALL interp_2d(x2, y2, f2, x1, y1, f1, dmx, dmy)

      DEALLOCATE(y2, x2, f2, dmx, dmy)
      RETURN
      END SUBROUTINE filter_2d
!----------------------------------------------------------------------
      SUBROUTINE smooth_dem(ft)
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) :: ft(:,:)
      REAL*8, ALLOCATABLE :: gt(:,:)
      INTEGER :: n1, n2, i, j

      n1 = SIZE(ft,1)
      n2 = SIZE(ft,2)
      ALLOCATE( gt(n1,n2) )
      gt = ft

      j = 1
      DO i = 2, n1-1
        gt(i,j) = 0.25D0 * ( ft(i-1,j) + 2.D0 * ft(i,j) + ft(i+1,j) )
      END DO
      j = n2
      DO i = 2, n1-1
        gt(i,j) = 0.25D0 * ( ft(i-1,j) + 2.D0 * ft(i,j) + ft(i+1,j) )
      END DO
      DO j = 2, n2-1
        i = 1
        gt(i,j) = 0.25D0 * ( ft(i,j-1) + 2.D0 * ft(i,j) + ft(i,j+1) )
        DO i = 2, n1-1
          gt(i,j) = 0.25D0 * ( ft(i,j-1) + 2.D0 * ft(i,j) + ft(i,j+1) )
        END DO
        i = n1
        gt(i,j) = 0.25D0 * ( ft(i,j-1) + 2.D0 * ft(i,j) + ft(i,j+1) )
      END DO

      ft = gt
      
      DEALLOCATE( gt )
      RETURN
      END SUBROUTINE smooth_dem
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
      SUBROUTINE set_profile
!
      USE control_flags, ONLY: job_type, lpr
      USE grid, ONLY: x, xb, y, yb, z, zb, iv, jv, kv

      IMPLICIT NONE
      INTEGER :: i,j,k,ijk
      LOGICAL, ALLOCATABLE :: dummy(:)
      REAL*8, ALLOCATABLE  :: topo(:)
      REAL*8, ALLOCATABLE  :: topo2d(:,:)
      REAL*8 :: transl_z = 0.D0
!
      ALLOCATE (dummy(ntot))
!
! ... Allocate and initialize 'dist' (defines implicitly the profile).
!
      ALLOCATE (dist(ntot))
      dist = 1.0D10
!
      IF (job_type == '2D') THEN
        !
        ALLOCATE(topo(nx))
        ALLOCATE(next(nx))
        ALLOCATE(ord(nx))

        ! ... Interpolate the topography on the cell top 
        ! ... (centered horizontally)
        !
        CALL interpolate_profile(x, zb, topo, dummy)
        
        ! ... Translate vertically the mesh to minimize the
        ! ... number of topographic cells
        !
        transl_z = MINVAL(topo) 
        IF( mpime == root ) THEN
          WRITE(6,*) 'Translating mesh vertically'
          WRITE(6,*) 'Minimum topographic quota: ', transl_z
        END IF
        z  = z  + transl_z
        zb = zb + transl_z
        
        ! ... Re-set the 'ord' array and set the implicit profile
        !
        IF (lpr > 1 .AND. mpime == root) &
          WRITE(7,*) 'topographic ordinates'
        DO i = 1, nx
          DO k = 1, nz
            IF (zb(k) <= topo(i)) ord(i) = k  
            ijk = i + (k-1) * nx
            dist(ijk) = z(k) - topo(i)
          END DO
          IF (lpr > 1 .AND. mpime == root) &
            WRITE(7,*) i, ord(i)
        END DO
        !
        DEALLOCATE(topo)
        !
      ELSE IF (job_type == '3D') THEN
        !
        ALLOCATE(topo2d(nx,ny))
        ALLOCATE(nextx(nx))
        ALLOCATE(nexty(ny))
        ALLOCATE(ord2d(nx,ny))

        CALL interpolate_dem(x, y, zb, topo2d, dummy)
        !
        ! ... Translate vertically the numerical mesh to minimize the
        ! ... number of topographic cells
        !
        transl_z = MINVAL(topo2d)
        IF( mpime == root ) THEN
          WRITE(6,*) 'Translating mesh vertically'
          WRITE(6,*) 'Minimum topographic quota: ', transl_z
        END IF
        z  = z  + transl_z
        zb = zb + transl_z
        
        ! ... Reset the 'ord2d' array and set the implicit profile
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
      CALL set_flag3

      DEALLOCATE (dummy)

      RETURN
      END SUBROUTINE set_profile
!----------------------------------------------------------------------
      SUBROUTINE interpolate_profile(cx, cz, topo, ff)
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
      CALL interp_1d(xtop,ztop,cx,topo,next)

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
      IMPLICIT NONE

      REAL*8, INTENT(IN), DIMENSION(:) :: cx, cy, cz
      REAL*8, INTENT(OUT), DIMENSION(:,:) :: topo2d
      LOGICAL, INTENT(OUT), DIMENSION(:) :: ff

      INTEGER i,j,k,ijk
      INTEGER itopo

      ff = .FALSE.

      CALL interp_2d(xtop, ytop, ztop2d, cx, cy, topo2d, nextx, nexty)
 
      ! ... Topography is expressed in centimeters
      !
      DO j=1,ny
        DO i=1,nx
          topo2d(i,j) = topo2d(i,j) * 1.D2
          itopo = NINT(topo2d(i,j))
          topo2d(i,j) = itopo * 1.D-2
        END DO
      END DO

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
      SUBROUTINE interp_1d(x1, f1, x2, f2, t)
      IMPLICIT NONE

      REAL*8, INTENT(IN), DIMENSION(:) :: x1, x2, f1
      REAL*8, INTENT(OUT), DIMENSION(:) :: f2
      INTEGER, INTENT(OUT), DIMENSION(:) :: t
      INTEGER :: i, k, l, n, n1x, n2x
      REAL*8 :: grad

      n1x = SIZE(x1)
      n2x = SIZE(x2)
!
! ... locate the grid points near the topographic points
! ... and interpolate linearly the profile  
!
      l = 1
      DO n = 1, n1x
        DO i = l, n2x

          IF (x1(n) >= x2(i)) THEN

            ! ... t(i) indicates the progressive number of the
            ! ... topographic point laying on the right of a grid center 'i'
            ! ... 'l' counts the grid points
            !
            t(i) = n
 
            IF (n == 1) THEN
              f2(i) = f1(1)
            ELSE
              grad = (f1(n)-f1(n-1))/(x1(n)-x1(n-1))
              f2(i) = f1(n-1) + (x2(i)-x1(n-1)) * grad
            END IF

            l=l+1
          ENDIF

        ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE interp_1d
!----------------------------------------------------------------------
      SUBROUTINE interp_2d(x1, y1, f1, x2, y2, f2, tx, ty)
      IMPLICIT NONE

      REAL*8, INTENT(IN), DIMENSION(:) :: x1, y1, x2, y2
      REAL*8, INTENT(IN), DIMENSION(:,:) :: f1
      REAL*8, INTENT(OUT), DIMENSION(:,:) :: f2
      INTEGER, INTENT(OUT), DIMENSION(:) :: tx, ty

      REAL*8 :: dist1y,dist2y,dist1x,dist2x,alpha,beta
      INTEGER :: i,j,k,h,l,ii,jj
      INTEGER :: n1x, n2x, n1y, n2y
      INTEGER :: tp1, tp2

      n1x = SIZE(x1)
      n1y = SIZE(y1)
      n2x = SIZE(x2)
      n2y = SIZE(y2)
        
!C============================================
!C===    trova le posizioni dei nodi      ====
!C===    della nuova griglia rispetto     ====
!C===    alla griglia iniziale            ====
!C============================================

      l=1
      DO i = 1, n1x
        DO ii = l, n2x
      
!============================================
!===    cerca i nodi della griglia che    ===
!===    che stanno a sx. di xtop(i)       ===
!============================================
    
          IF (x1(i) >= x2(ii)) THEN
            tx(ii)=i
            l=l+1
          ENDIF
        ENDDO
      ENDDO

      l=1
      DO j = 1, n1y
        DO jj = l, n2y

!C============================================
!C===    cerca i nodi della griglia che    ===
!C===    che stanno a sotto  ytop(i)       ===
!C============================================

          IF (y1(j) >= y2(jj)) THEN
            ty(jj)=j
            l=l+1
          ENDIF
        ENDDO
      ENDDO


! il nodo della nuova griglia di indici (i,j) sara' allora
! contenuto nel rettangolo con i vertici con indici:
! P1=tx(i),ty(j)
! P2=tx(i-1),ty(j)
! P3=tx(i-1),ty(j-1)
! P4=tx(i),ty(j-1)

! sulla nuova griglia interpoliamo le quote di input ztop per 
! ottenere la quota coorZ nel punto di indici (i,j)
 

! interpolazione bilineare sui nodi interni (1<i<nodiGRIDx, 1<j<nodigGRIDy)
! utilizzando le quote nei punti P1,..,P4 definiti sopra

        DO i=1,n2x

           dist1x = x2(i) - x1(tx(i)-1)
           dist2x = x1(tx(i)) - x2(i)
           alpha  = dist1x/(dist1x+dist2x)

           DO j=1,n2y
           
              dist1y = y2(j) - y1(ty(j)-1)
              dist2y = y1(ty(j)) - y2(j)
              beta   = dist1y/(dist1y+dist2y)

              tp1    = alpha * f1(tx(i),ty(j))   + &
                       (1.D0 - alpha) * f1(tx(i)-1,ty(j))
              tp2    = alpha * f1(tx(i),ty(j)-1) + &
                       (1.D0 - alpha) * f1(tx(i)-1,ty(j)-1)
              f2(i,j) = beta * tp1 + (1.D0 - beta) * tp2
 
           ENDDO
        ENDDO
!      
      RETURN
      END SUBROUTINE interp_2d
!----------------------------------------------------------------------
      SUBROUTINE write_profile

      USE control_flags, ONLY: job_type
      USE grid, ONLY: x, xb, y, yb, z, zb
      USE grid, ONLY: iob, fl, bottom
!
      IMPLICIT NONE
      INTEGER :: n, i, j, k, ijk
!
! ... Write the georeferenced mesh coordinates
!
      OPEN(17,FILE='mesh.dat')
      WRITE(17,*) 'Georeferenced x-y mesh'

      WRITE(17,*) 'x'
      WRITE(17,17) x
      WRITE(17,*) 'xb'
      WRITE(17,17) xb
      WRITE(17,*) 'y'
      WRITE(17,17) y
      WRITE(17,*) 'yb'
      WRITE(17,17) yb
      WRITE(17,*) 'z'
      WRITE(17,17) z
      WRITE(17,*) 'zb'
      WRITE(17,17) zb
      CLOSE(17)

 17   FORMAT(5(F20.6))
!
! ... Write out the implicit profile
!
      IF (mpime == root) THEN
        OPEN(UNIT=14,FILE='improfile.dat',STATUS='UNKNOWN')
        WRITE(14,fmt='(F9.2)') dist
        CLOSE(14)
      END IF
!
! ... Control that the cells below the flag = 5 blocks
! ... have flag = 'bottom'
!
      IF( job_type == '2D') THEN
        DO n = 1, no
          IF (iob(n)%typ == 5) THEN
            DO i = iob(n)%xlo, iob(n)%xhi
              DO k = 1, iob(n)%zlo - 1
                ijk = i + (k-1) * nx
                fl(ijk) = bottom
              END DO
            END DO
          END IF
        END DO
      ELSE IF (job_type == '3D') THEN
        DO n = 1, no
          IF (iob(n)%typ == 5) THEN
            DO i = iob(n)%xlo, iob(n)%xhi
              DO j = iob(n)%ylo, iob(n)%yhi
                DO k = 1, iob(n)%zlo - 1
                  ijk = i + (j-1) * nx + (k-1) * nx * ny
                  fl(ijk) = bottom
                END DO
              END DO
            END DO
          END IF
        END DO
      END IF
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
! 
      RETURN
      END SUBROUTINE write_profile
!----------------------------------------------------------------------
      SUBROUTINE set_flag3
!
! ... Set cell-flag = 3 in cells laying below the topography
!
      USE control_flags, ONLY: lpr, job_type
      USE grid, ONLY: fl
      IMPLICIT NONE

      INTEGER :: i, j, k, ijk
      INTEGER :: q
!
      IF( job_type == '2D') THEN
        DO i = 1, nx
          q = ord(i)
          DO k = 1, nz
            ijk = i + (k-1) * nx
            IF (q >= k) THEN
                    IF( fl(ijk)/=5 ) fl(ijk) = 3
            ELSE
              EXIT
            END IF
          END DO
        END DO
      ELSE IF( job_type == '3D') THEN
        DO j = 1, ny
          DO i = 1, nx
            q = ord2d(i,j)
            DO k = 1, nz
              ijk = i + (j-1) * nx + (k-1) * nx * ny
              IF (q >= k) THEN
                      IF( fl(ijk)/=5 ) fl(ijk) = 3
              ELSE
                EXIT
              END IF
            END DO
          END DO
        END DO
      END IF
!
      RETURN
      END SUBROUTINE set_flag3
!----------------------------------------------------------------------
      END MODULE volcano_topography
!----------------------------------------------------------------------
