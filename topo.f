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
      INTEGER :: noditop, noditopx, noditopy
!
      INTEGER :: itp, iavv
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
        CALL read_dem_ascii
        IF (iavv >= 1) CALL averaged_volcano
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
      USE grid, ONLY: topography
      USE parallel, ONLY: mpime, root
      IMPLICIT NONE

      INTEGER :: n

      CHARACTER(LEN=80) :: topo_file
!
! ... Processor 'root' reads the topography defined on a 
! ... (generic) rectilinear mesh "xtop". Then broadcast
! ... to all processors.
!
      topo_file = TRIM(topography)
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
      USE grid, ONLY: topography
      USE parallel, ONLY: mpime, root
      IMPLICIT NONE

      CHARACTER(LEN=80) :: topo_file
      REAL    :: xll, yll
      REAL    :: dd
      INTEGER :: noval
      INTEGER :: elevation
      INTEGER :: p
      INTEGER :: i, j
!
! ... Processor 'root' reads the standard header of the
! ... ASCII DEM format, then broadcasts it to all processors.
!
      IF (mpime == root) THEN
        topo_file = TRIM(topography)
        OPEN(UNIT=3, FILE=topo_file, STATUS='OLD')
        WRITE(6,*) 'Reading topography file: ', topo_file
        READ(3,*) noditopx
        READ(3,*) noditopy
        READ(3,*) xll
        READ(3,*) yll
        READ(3,*) dd
        READ(3,*) noval
      END IF
!
      CALL bcast_integer(noditopx,1,root)
      CALL bcast_integer(noditopy,1,root)
      CALL bcast_real(xll,1,root)
      CALL bcast_real(yll,1,root)
      CALL bcast_real(dd,1,root)
      CALL bcast_integer(noval,1,root)
!
      vdem%nx           = noditopx
      vdem%ny           = noditopy
      vdem%xcorner      = xll
      vdem%ycorner      = yll
      vdem%cellsize     = dd
      vdem%nodata_value = noval
!
      ALLOCATE(ztop2d(vdem%nx,vdem%ny))
      ALLOCATE(xtop(vdem%nx))
      ALLOCATE(ytop(vdem%ny))
!
! ... Compute the UTM coordinates of each element of the topography
!
      DO i = 1, vdem%nx
        xtop(i) = vdem%xcorner + (i-1) * vdem%cellsize
      END DO
      DO j = vdem%ny, 1, -1
        ytop(j) = vdem%ycorner - (vdem%ny - j) * vdem%cellsize
      END DO
!
! ... Processor 'root' reads the matrix of the elevation
! ... ( in centimeters ) and broadcasts.
!
      IF (mpime == root) THEN
        DO j = vdem%ny, 1, -1
          DO i = 1, vdem%nx
            READ(3,*) elevation
            ztop2d(i,j) = DBLE(elevation) / 100.D0
          END DO
        END DO
      END IF
!
      CALL bcast_real(ztop2d,noditopx*noditopy,root)
!
      RETURN
      END SUBROUTINE read_dem_ascii
!----------------------------------------------------------------------
      SUBROUTINE averaged_volcano
      USE grid, ONLY: center_x, center_y
      IMPLICIT NONE
      INTEGER :: i, j, icenter, jcenter, distance, m
      INTEGER :: dms
      REAL*8, ALLOCATABLE :: av_quota(:)
      INTEGER, ALLOCATABLE :: nk(:)
!
      dms = vdem%nx**2 + vdem%ny**2
      ALLOCATE( av_quota( 0:dms ) )
      ALLOCATE( nk( 0:dms ) )
      av_quota = 0.D0
      nk       = 0
!
      ! ... Find the topographic element much closer
      ! ... to the vent center
      !
      DO i = 1, vdem%nx
        IF (xtop(i) <= center_x) icenter = i
      END DO
      DO j = 1, vdem%ny
        IF (ytop(j) <= center_y) jcenter = j
      END DO
!      
      ! ... The vent center must coincide with the
      ! ... topographic element
      !
      center_x = xtop(icenter)
      center_y = ytop(jcenter)
!
      ! ... Average the topography over axisymmetrix layers
      !
      DO j = 1, vdem%ny
        DO i = 1, vdem%nx
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
      ! ... Assign the new averaged values to the
      ! ... topographic arrays
      !
      DO j = 1, vdem%ny
        DO i = 1, vdem%nx
          distance = (i - icenter)**2 + (j - jcenter)**2
          ztop2d(i,j) = av_quota(distance)
        END DO
      END DO
!
      DEALLOCATE(av_quota)
      DEALLOCATE(nk)
!
      RETURN
      END SUBROUTINE averaged_volcano
!----------------------------------------------------------------------
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
! ... Translate all mesh horizontally
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
        CALL interpolate_2d(x, zb, topo, dummy)
        
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
      SUBROUTINE interpolate_2d(cx, cz, topo, ff)
!
! ... interpolate the topographic profile on a given
! ... computational mesh (either centered or staggered)
! ... as defined by the arrays (cx, cy, cz)

      IMPLICIT NONE

      REAL*8, INTENT(IN), DIMENSION(:) :: cx, cz
      LOGICAL, INTENT(OUT), DIMENSION(:) :: ff
      REAL*8, INTENT(OUT), DIMENSION(:) :: topo

      INTEGER :: i, k, ijk, l, n
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

            ! ... Topography is expressed in centimeters
            !
            topo(i) = topo(i) * 1.D2
            itopo = NINT(topo(i))
            topo(i) = itopo * 1.D-2

            l=l+1
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
      END SUBROUTINE interpolate_2d
!----------------------------------------------------------------------
      SUBROUTINE interpolate_dem(cx, cy, cz, topo2d, ff)
      USE dimensions, ONLY: nx, ny, nz
      IMPLICIT NONE

      REAL*8, INTENT(IN), DIMENSION(:) :: cx, cy, cz
      REAL*8, INTENT(OUT), DIMENSION(:,:) :: topo2d
      LOGICAL, INTENT(OUT), DIMENSION(:) :: ff

      REAL*8 xmin,xmax,ymin,ymax,zmin,zmax !range dei valori della top.
      REAL*8 ratio

      REAL*8 dist1y,dist2y,dist1x,dist2x,alpha,beta
      INTEGER i,j,k,h,l,ii,jj
      INTEGER ijk, tp1, tp2
      INTEGER itopo

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
 
              ! ... Topography is expressed in centimeters
              !
              topo2d(i,j) = topo2d(i,j) * 1.D2
              itopo = NINT(topo2d(i,j))
              topo2d(i,j) = itopo * 1.D-2

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
      INTEGER :: quota
!
      IF( job_type == '2D') THEN
        DO i = 1, nx
          quota = ord(i)
          DO k = 1, nz
            ijk = i + (k-1) * nx
            IF (quota >= k) THEN
                    IF( fl(ijk)/=5 ) fl(ijk) = 3
            ELSE
              EXIT
            END IF
          END DO
        END DO
      ELSE IF( job_type == '3D') THEN
        DO j = 1, ny
          DO i = 1, nx
            quota = ord2d(i,j)
            DO k = 1, nz
              ijk = i + (j-1) * nx + (k-1) * nx * ny
              IF (quota >= k) THEN
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
