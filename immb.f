!----------------------------------------------------------------------
      MODULE immersed_boundaries
!----------------------------------------------------------------------
      USE dimensions, ONLY: nx, ny, nz, ntot, ntr
      USE parallel, ONLY: mpime, root
      USE volcano_topography, ONLY: xtop, ytop, ztop, ztop2d

      IMPLICIT NONE
!
      TYPE point
        REAL*8 :: x
        REAL*8 :: y
        REAL*8 :: z
      END TYPE point
!
! ... Parameters identifying a forcing point: (i,j,k) are the 
! ... (x,y,z) discrete coordinates, (int) is the type of interpolation
! ... to be done:
! ... In 2D: -10= extrap z -3=lin sx/top, -2=linsx, -1=bilsx, 
! ...        0=lintop, 1=bildx, 2=lindx, 3=lin dx/top, 10=extrap x), 
! ... In 3D: the interpolation is defined by the number of neighbours
! ...        and identified by the binary numbers form 0 to 255
! ... (nsl) are the coordinates of the noslip point
!
      TYPE forcing_point
        INTEGER :: i
        INTEGER :: j
        INTEGER :: k
        INTEGER :: int
        TYPE(point) :: nsl
        REAL*8  :: vel
      END TYPE forcing_point
!
      TYPE(forcing_point), ALLOCATABLE :: fptx(:), fpty(:), fptz(:)
      TYPE(forcing_point), ALLOCATABLE :: forx(:), fory(:), forz(:)
!
! ... Conditional array to be used in loops:
! ... this value is .TRUE. when force has to be computed at a given location
!
      LOGICAL, ALLOCATABLE :: forcex(:)
      LOGICAL, ALLOCATABLE :: forcey(:)
      LOGICAL, ALLOCATABLE :: forcez(:)
      LOGICAL, ALLOCATABLE :: forced(:)
!
! ... Topographic elevation at the mesh points
! ... (centered and staggered)
!
      REAL*8, ALLOCATABLE  :: topo_c(:)
      REAL*8, ALLOCATABLE  :: topo_x(:)
      REAL*8, ALLOCATABLE  :: topo2d_c(:,:)
      REAL*8, ALLOCATABLE  :: topo2d_x(:,:)
      REAL*8, ALLOCATABLE  :: topo2d_y(:,:)
!
! ... This function gives the increasing number of local forcing points 
      INTEGER, ALLOCATABLE :: numx(:), numy(:), numz(:)
      INTEGER :: nfp
!           
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE set_forcing
      USE control_flags, ONLY: job_type, lpr, immb, itp
      USE volcano_topography, ONLY: grid_locations, vertical_shift
      USE volcano_topography, ONLY: interpolate_2d, interpolate_dem
      USE volcano_topography, ONLY: ord, next
      USE volcano_topography, ONLY: ord2d, nextx, nexty
      USE volcano_topography, ONLY: cx, cy, cz
      USE grid, ONLY: z

      IMPLICIT NONE
      INTEGER :: p, nfpx, nfpy, nfpz
      INTEGER :: i,j,k,ijk

      LOGICAL, ALLOCATABLE :: dummy(:)
      ALLOCATE (dummy(ntot))
!
! ... If Immersed Boundaries are used, identify the forcing points
! ... and set interpolation parameters
!
      !
      ! ... Allocate the logical arrays that are used to 
      ! ... identify the forcing points
      !
      ALLOCATE(forcex(ntot))
      IF (job_type == '3D') ALLOCATE(forcey(ntot))
      ALLOCATE(forcez(ntot))

      IF (job_type == '2D') THEN

        ALLOCATE(topo_c(nx))
        ALLOCATE(topo_x(nx))
        !
        ! ... interpolate the topography on x-staggered mesh
        !
        CALL grid_locations(1,0,0)
        CALL interpolate_2d(topo_x,forcex)
        nfpx = COUNT(forcex)
        ALLOCATE(fptx(nfpx))
        !
        ! ... Forcing along x
        CALL forcing2d(fptx,topo_x)
        !
        ! ... interpolate the topography on z-staggered mesh
        !
        CALL grid_locations(0,0,1)
        CALL interpolate_2d(topo_c,forcez)
        nfpz = COUNT(forcez)
        ALLOCATE(fptz(nfpz))
        !
        ! ... Forcing along z
        CALL forcing2d(fptz,topo_c)
        !
        ! ... Merge forcing points subsets
        !
        nfp = COUNT( forcex .OR. forcez )
        ALLOCATE (forx (nfp) )
        ALLOCATE (forz (nfp) )

        CALL grid_locations(1,0,0)
        CALL interpolate_2d(topo_x, dummy)
        forx(1:nfpx) = fptx(1:nfpx)
        CALL extrafx2d(forx, nfpx, topo_x)

        CALL grid_locations(0,0,1)
        CALL interpolate_2d(topo_c, dummy)
        forz(1:nfpz) = fptz(1:nfpz)
        CALL extrafz2d(forz, nfpz, topo_c)
        
        DEALLOCATE(fptx)
        DEALLOCATE(fptz)

      ELSE IF (job_type == '3D') THEN

        ALLOCATE(topo2d_c(nx,ny))
        ALLOCATE(topo2d_x(nx,ny))
        ALLOCATE(topo2d_y(nx,ny))
        !
        ! ... interpolate the topography on x-staggered mesh
        !
        CALL grid_locations(1,0,0)
        CALL interpolate_dem(topo2d_x, forcex)
        nfpx = COUNT(forcex)
        ALLOCATE(fptx(nfpx))
        !
        ! ... Forcing along x
        CALL forcing3d(fptx,topo2d_x)
        !
        ! ... interpolate the topography on y-staggered mesh
        !
        CALL grid_locations(0,1,0)
        CALL interpolate_dem(topo2d_y, forcey)
        nfpy = COUNT(forcey)
        ALLOCATE(fpty(nfpy))
        !
        ! ... Forcing along y
        CALL forcing3d(fpty,topo2d_y)
        !
        ! ... interpolate the topography on z-staggered mesh
        !
        CALL grid_locations(0,0,1)
        CALL interpolate_dem(topo2d_c, forcez)
        nfpz = COUNT(forcez)
        ALLOCATE(fptz(nfpz))
        !
        ! ... Forcing along z
        CALL forcing3d(fptz,topo2d_c)
        !
        ! ... Merge forcing points subsets
        !
        nfp = COUNT( forcex .OR. forcey .OR. forcez )
        ALLOCATE (forx (nfp) )
        ALLOCATE (fory (nfp) )
        ALLOCATE (forz (nfp) )

        CALL grid_locations(1,0,0)
        CALL interpolate_dem(topo2d_x, dummy)
        forx(1:nfpx) = fptx(1:nfpx)
        CALL extrafx3d(forx, nfpx, topo2d_x)

        CALL grid_locations(0,1,0)
        CALL interpolate_dem(topo2d_y, dummy)
        fory(1:nfpy) = fpty(1:nfpy)
        CALL extrafy3d(fory, nfpy, topo2d_y)

        CALL grid_locations(0,0,1)
        CALL interpolate_dem(topo2d_c, dummy)
        forz(1:nfpz) = fptz(1:nfpz)
        CALL extrafz3d(forz, nfpz, topo2d_c)
!
        DEALLOCATE(fptx)
        DEALLOCATE(fpty)
        DEALLOCATE(fptz)

      END IF

      IF (mpime == root) THEN
        OPEN(UNIT=15,FILE='forx.dat',STATUS='UNKNOWN')
        IF (job_type == '3D') &
          OPEN(UNIT=16,FILE='fory.dat',STATUS='UNKNOWN')
        OPEN(UNIT=17,FILE='forz.dat',STATUS='UNKNOWN')
        DO p = 1, nfp
          WRITE(15,33) p, forx(p)
          IF (job_type == '3D') WRITE(16,33) p, fory(p)
          WRITE(17,33) p, forz(p)
        END DO
        CLOSE(15)
        CLOSE(16)
        CLOSE(17)
      END IF
 33   FORMAT(5(I6),4(F18.3))

      IF (job_type == '2D') THEN
        DEALLOCATE (topo_c, topo_x)
      ELSE IF (job_type == '3D') THEN
        DEALLOCATE (topo2d_c, topo2d_x, topo2d_y)
      END IF
      DEALLOCATE (dummy)
!
      RETURN
      END SUBROUTINE set_forcing
!----------------------------------------------------------------------
      SUBROUTINE forcing2d(fpt, topo)
      USE volcano_topography, ONLY: cx, cy, cz
      USE volcano_topography, ONLY: ord, next
!
! ... locate the forcing points and the type of interpolation
! ... locate the points where no-slip condition is applied

      IMPLICIT NONE
      INTEGER :: i, fp
      REAL*8, DIMENSION(:), INTENT(IN) :: topo
      TYPE(forcing_point), DIMENSION(:), INTENT(OUT) :: fpt
      REAL*8 :: t,s
!
! ... on the -- first grid point -- interpolate   
! ... linearly along the vertical direction  
!
      i = 1
      fp = 1

      fpt(fp)%int = 0	
      fpt(fp)%nsl%x = cx(i)
      fpt(fp)%nsl%z = topo(i)
      fpt(fp)%i  = i
      fpt(fp)%k  = ord(i)
!
! ... loop over grid points in x-direction
!
      DO i=2,nx-1

        IF ((ord(i)-ord(i-1) > 0).AND.(ord(i+1)-ord(i) >= 0)) THEN
          
          CALL increasing_profile(i,fp,fpt)

        ELSEIF (ord(i+1)-ord(i) < 0) THEN

          CALL decreasing_profile(i,fp,fpt)

        ELSE 

          CALL stationnary_profile(i,fp,fpt)

        END IF

      END DO
!
! ... on the -- last grid point -- interpolate    
! ... linearly along the vertical direction  
!
      i = nx
      fp = fp + 1

      fpt(fp)%int = 0	
      fpt(fp)%nsl%x = cx(i)
      fpt(fp)%nsl%z = topo(i)
      fpt(fp)%i  = i
      fpt(fp)%k  = ord(i)

      RETURN
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE increasing_profile(i,fp,fpt)
      USE volcano_topography, ONLY: cx, cy, cz

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(INOUT) :: fp
      TYPE(forcing_point), DIMENSION(:), INTENT(INOUT) :: fpt

      REAL*8 :: norm(2)
      REAL*8 :: grad, s
      REAL*8 :: dxt, dzt, dxs, dzs
      INTEGER :: n, k, err
!
! ... diagonal interpolation 
!
      fp = fp + 1

      fpt(fp)%int = -3
      fpt(fp)%i  = i
      fpt(fp)%k  = ord(i-1)
      
      norm(1) = cx(i-1) - cx(i)
      norm(2) = cz(ord(i-1)+1) - cz(ord(i-1))
!
! ... find the no-slip point on the profile 
!
      DO n=next(i-1),next(i)
	
        dxt = xtop(n) - xtop(n-1)
        dzt = ztop(n) - ztop(n-1)
        dxs = xtop(n) - cx(i)
        dzs = ztop(n) - cz(ord(i-1))

        s = sol2(norm(1), dxt, norm(2), dzt, dxs, dzs, err)
        IF (err > 0) WRITE(6,*) 'Error in fp: ', fp

        IF ((s >= 0).AND.(s <= 1)) THEN
          fpt(fp)%nsl%x = s*xtop(n-1) + (1-s)*xtop(n)
          fpt(fp)%nsl%z = s*ztop(n-1) + (1-s)*ztop(n)
        ENDIF
      ENDDO
!
! ... loop over grid points with same x-coordinate
! ... (left linear interpolation)
!
      DO k = ord(i-1) + 1 , ord(i) - 1
        fp = fp + 1
        fpt(fp)%int = -2
        fpt(fp)%i  = i
        fpt(fp)%k  = k
!
! ... find the no-slip point on the profile
!
        DO n=next(i-1),next(i)
          IF ((ztop(n-1) <= cz(k)).AND. (ztop(n) >= cz(k))) THEN
            grad = (xtop(n)-xtop(n-1))/(ztop(n)-ztop(n-1))
            fpt(fp)%nsl%x = xtop(n-1) + (cz(k)-ztop(n-1)) * grad
            fpt(fp)%nsl%z = cz(k)
          ENDIF
        ENDDO
          
      ENDDO
!
! ... bilinear interpolation
! ... on the last (upper) point on the same x-coord
!
      fp = fp + 1
      fpt(fp)%int = -1
      fpt(fp)%i  = i
      fpt(fp)%k  = ord(i)

      norm(1) = topo(i-1)-topo(i)
      norm(2) = cx(i)-cx(i-1)
!
! ... find the no-slip point on the profile
!
                    
      DO n=next(i-1),next(i)

        dxt = xtop(n)-xtop(n-1)
        dzt = ztop(n)-ztop(n-1)
        dxs = xtop(n) - cx(i)
        dzs = ztop(n) - cz(ord(i))

        s = sol2( norm(1), dxt, norm(2), dzt, dxs, dzs , err)
        IF (err > 0) WRITE(6,*) 'Error in fp: ', fp

        IF ((s >= 0).AND.(s <= 1))	THEN
          fpt(fp)%nsl%x = s*xtop(n-1) + (1-s)*xtop(n)
          fpt(fp)%nsl%z = s*ztop(n-1) + (1-s)*ztop(n)
        ENDIF
      ENDDO
      
      END SUBROUTINE increasing_profile
!----------------------------------------------------------------------
      SUBROUTINE decreasing_profile(i,fp,fpt)
      USE volcano_topography, ONLY: cx, cy, cz

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(INOUT) :: fp
      TYPE(forcing_point), DIMENSION(:), INTENT(INOUT) :: fpt

      REAL*8 :: norm(2)
      REAL*8 :: grad, s
      REAL*8 :: dxt, dzt, dxs, dzs
      INTEGER :: n, k, err
!
! ... bilinear interpolation
! ... on the first (upper) point
!
      fp = fp + 1

      fpt(fp)%int=1
      fpt(fp)%i  = i
      fpt(fp)%k  = ord(i)

      norm(1) = topo(i) - topo(i+1)
      norm(2) = cx(i+1) - cx(i)
!
! ... find the no-slip point on the profile
!
      DO n=next(i),next(i+1)

        dxt = xtop(n) - xtop(n-1)
        dzt = ztop(n) - ztop(n-1)
        dxs = xtop(n) - cx(i)
        dzs = ztop(n) - cz(ord(i))

	s = sol2( norm(1), dxt, norm(2), dzt, dxs, dzs, err )
        IF (err > 0) WRITE(6,*) 'Error in fp: ', fp

        IF ((s >= 0).AND.(s <= 1))	THEN
          fpt(fp)%nsl%x = s*xtop(n-1) + (1-s)*xtop(n)
          fpt(fp)%nsl%z = s*ztop(n-1) + (1-s)*ztop(n)
        ENDIF
      ENDDO
!
! ... loop over grid points with same x-coordinate
! ... (right linear interpolation)
!
      DO k=ord(i)-1,ord(i+1)+1,-1

       fp = fp + 1

       fpt(fp)%int = 2
       fpt(fp)%i  = i
       fpt(fp)%k  = k
!
! ... find the no-slip point on the profile
!
        DO n=next(i),next(i+1)
          IF ((ztop(n-1) >= cz(k)).AND.(ztop(n) <= cz(k))) THEN
            grad = (xtop(n)-xtop(n-1))/(ztop(n)-ztop(n-1))
	    fpt(fp)%nsl%x = xtop(n-1) + (cz(k)-ztop(n-1)) * grad
	    fpt(fp)%nsl%z = cz(k)
	  ENDIF
	ENDDO

      ENDDO
!
! ... diagonal interpolation 
!
      fp = fp + 1

      fpt(fp)%int = 3
      fpt(fp)%i  = i
      fpt(fp)%k  = ord(i+1)

      norm(1)=cx(i+1)-cx(i)
      norm(2)=cz(ord(i+1)+1)-cz(ord(i+1))
!
! ... find the no-slip point on the profile
!
      DO n=next(i),next(i+1)
 
        dxt = xtop(n) - xtop(n-1)
        dzt = ztop(n) - ztop(n-1)
        dxs = xtop(n) - cx(i)
        dzs = ztop(n) - cz(ord(i+1))

	s = sol2( norm(1), dxt, norm(2), dzt, dxs, dzs, err )
        IF (err > 0) WRITE(6,*) 'Error in fp: ', fp

 	IF ((s >= 0).AND.(s <= 1))	THEN
	  fpt(fp)%nsl%x = s*xtop(n-1) + (1-s)*xtop(n)
	  fpt(fp)%nsl%z = s*ztop(n-1) + (1-s)*ztop(n)
 	ENDIF
      ENDDO
      
      END SUBROUTINE decreasing_profile
!----------------------------------------------------------------------
      SUBROUTINE stationnary_profile(i,fp,fpt)
      USE volcano_topography, ONLY: cx, cy, cz

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(INOUT) :: fp
      TYPE(forcing_point), DIMENSION(:), INTENT(INOUT) :: fpt

      fp = fp + 1

      fpt(fp)%int=0	

      fpt(fp)%i  = i
      fpt(fp)%k  = ord(i)

      fpt(fp)%nsl%x = cx(i)
      fpt(fp)%nsl%z = topo(i)

      END SUBROUTINE stationnary_profile
!----------------------------------------------------------------------
      FUNCTION sol1(a11,a12,a21,a22,b1,b2,info)

      REAL*8 a11,a12,a21,a22	! coefficienti della matrice
      REAL*8 b1,b2		! vettore dei termini noti
      REAL*8 sol1
      INTEGER, INTENT(OUT) :: info

      info = 0
      IF (a11*a22-a12*a21 == 0) info = 1
      sol1=(b1*a22-a12*b2)/(a11*a22-a12*a21)

      END FUNCTION sol1
!----------------------------------------------------------------------
      FUNCTION sol2(a11,a12,a21,a22,b1,b2, info)

      REAL*8 a11,a12,a21,a22	! coefficienti della matrice
      REAL*8 b1,b2		! vettore dei termini noti
      REAL*8 sol2
      INTEGER, INTENT(OUT) :: info
 
      info = 0
      IF (a11*a22-a12*a21 == 0) info = 1
      sol2=(a11*b2-b1*a21)/(a11*a22-a12*a21)

      END FUNCTION sol2
!----------------------------------------------------------------------
!
      END SUBROUTINE forcing2d
!
!----------------------------------------------------------------------
      SUBROUTINE extrafx2d(fx, nfpx, topo)
      USE volcano_topography, ONLY: cx, cy, cz, next
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nfpx
      TYPE(forcing_point), DIMENSION(:), INTENT(INOUT) :: fx
      REAL*8, DIMENSION(:), INTENT(IN) :: topo
      INTEGER :: ijk, i, j, k
      INTEGER :: npx, n
      REAL*8 :: grad
        
      npx = nfpx

      DO ijk = 1, ntot
!
! ... If the vertical velocity is forced but the horizontal
! ... velocity is not, impose the forcing horizontally with
! ... extrapolation ('int'= -10)
!
        IF (forcez(ijk).AND.(.NOT.forcex(ijk))) THEN

          k = ( ijk - 1 ) / nx + 1
          i = MOD( ( ijk - 1 ), nx) + 1

          npx = npx + 1
          fx(npx)%i = i
          fx(npx)%k = k
          !
          ! ... find the no-slip point on the profile
          DO n=next(i),next(i+1)
            IF ((ztop(n-1) >= cz(k)).AND.(ztop(n) <= cz(k))) THEN
              grad = (xtop(n)-xtop(n-1))/(ztop(n)-ztop(n-1))
  	      fx(npx)%nsl%x = xtop(n-1) + (cz(k)-ztop(n-1)) * grad
  	      fx(npx)%nsl%z = cz(k)
  	    ENDIF
          ENDDO

          IF (cz(k) >= topo(i)) THEN
            ! ... extrapolation is done if the point lays
            ! ... above the topography ...
            !
            fx(npx)%int = -10
          ELSE
            ! ... otherwise do nothing!
            !
            fx(npx)%int = -100
            fx(npx)%nsl%x = 0.D0
            fx(npx)%nsl%z = 0.D0
          END IF

        END IF
      END DO
!
      IF (npx/=nfp) CALL error('extrafx','check nfp',1)
!
      RETURN
      END SUBROUTINE extrafx2d
!----------------------------------------------------------------------
      SUBROUTINE extrafz2d(fz, nfpz, topo)
      USE volcano_topography, ONLY: cx, cy, cz, next
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nfpz
      TYPE(forcing_point), DIMENSION(:), INTENT(INOUT) :: fz
      REAL*8, DIMENSION(:), INTENT(IN) :: topo
      INTEGER :: ijk, i, j, k
      INTEGER :: npz
        
      npz = nfpz

      DO ijk = 1, ntot
        k = ( ijk - 1 ) / nx + 1
        i = MOD( ( ijk - 1 ), nx) + 1
        IF (forcex(ijk).AND.(.NOT.forcez(ijk))) THEN
          npz = npz + 1
          fz(npz)%i = i
          fz(npz)%k = k
          IF (cz(k) >= topo(i)) THEN
            fz(npz)%int = 10
            fz(npz)%nsl%x = cx(i)
            fz(npz)%nsl%z = topo(i)
          ELSE
            fz(npz)%int = 100
          END IF
        END IF
      END DO
!
      IF (npz/=nfp) CALL error('extrafz','check nfp',1)

      RETURN
      END SUBROUTINE extrafz2d
!----------------------------------------------------------------------
      SUBROUTINE forcing3d(fpt, topo2d)
      USE volcano_topography, ONLY: ord2d, next
      USE volcano_topography, ONLY: cx, cy, cz
!
! ... locate the forcing points and the type of interpolation
! ... locate the points where no-slip condition is applied
!
      IMPLICIT NONE
      INTEGER :: i, j, k, fp
      REAL*8, DIMENSION(:,:), INTENT(IN) :: topo2d
      TYPE(forcing_point), DIMENSION(:), INTENT(OUT) :: fpt
      REAL*8 :: g(8)
      INTEGER :: delta_i(8), delta_j(8)
      INTEGER :: gint
!
      fp = 0
!
! ... LINEAR interpolation criteria ('fpt%int'):
! ... 0: top
! ... 1: west
! ... 2: south-west
! ... 3: south
! ... 4: south-east
! ... 5: east
! ... 6: north-east
! ... 7: north
! ... 8: north-west
!
      delta_i(1) = -1
      delta_j(1) = 0
      delta_i(2) = -1
      delta_j(2) = -1
      delta_i(3) = 0
      delta_j(3) = -1
      delta_i(4) = 1
      delta_j(4) = -1
      delta_i(5) = 1
      delta_j(5) = 0
      delta_i(6) = 1
      delta_j(6) = 1
      delta_i(7) = 0
      delta_j(7) = 1
      delta_i(8) = -1
      delta_j(8) = 1

      i = 1
      DO j = 1, ny
        k = ord2d(i,j)

        fp = fp + 1 
        fpt(fp)%i  = i
        fpt(fp)%j  = j
        fpt(fp)%k  = k
        fpt(fp)%int = 0
        fpt(fp)%nsl%x = cx(i)
        fpt(fp)%nsl%y = cy(j)
        fpt(fp)%nsl%z = topo2d(i,j)

      END DO

      i = nx
      DO j = 1, ny 
        k = ord2d(i,j)

        fp = fp + 1 
        fpt(fp)%i  = i
        fpt(fp)%j  = j
        fpt(fp)%k  = k
        fpt(fp)%int = 0
        fpt(fp)%nsl%x = cx(i)
        fpt(fp)%nsl%y = cy(j)
        fpt(fp)%nsl%z = topo2d(i,j)

      END DO

      DO i = 2, nx - 1

        j = 1
        k = ord2d(i,j)

        fp = fp + 1 
        fpt(fp)%i  = i
        fpt(fp)%j  = j
        fpt(fp)%k  = k
        fpt(fp)%int = 0
        fpt(fp)%nsl%x = cx(i)
        fpt(fp)%nsl%y = cy(j)
        fpt(fp)%nsl%z = topo2d(i,j)

        j = ny
        k = ord2d(i,j)

        fp = fp + 1 
        fpt(fp)%i  = i
        fpt(fp)%j  = j
        fpt(fp)%k  = k
        fpt(fp)%int = 0
        fpt(fp)%nsl%x = cx(i)
        fpt(fp)%nsl%y = cy(j)
        fpt(fp)%nsl%z = topo2d(i,j)

        DO j = 2, ny - 1
          DO k = 1, ord2d(i,j) - 1
            g = -1
            IF( ord2d(i-1,j)   < k ) g(1) = gamma(i,j,k, 1)
            IF( ord2d(i-1,j-1) < k ) g(2) = gamma(i,j,k, 2)
            IF( ord2d(i,j-1)   < k ) g(3) = gamma(i,j,k, 3)
            IF( ord2d(i+1,j-1) < k ) g(4) = gamma(i,j,k, 4)
            IF( ord2d(i+1,j)   < k ) g(5) = gamma(i,j,k, 5)
            IF( ord2d(i+1,j+1) < k ) g(6) = gamma(i,j,k, 6)
            IF( ord2d(i,j+1)   < k ) g(7) = gamma(i,j,k, 7)
            IF( ord2d(i-1,j+1) < k ) g(8) = gamma(i,j,k, 8)

            IF (MAXVAL(g) > -1) THEN

              gint = mingam(g) 

              fp = fp + 1 
              fpt(fp)%i  = i
              fpt(fp)%j  = j
              fpt(fp)%k  = k
              fpt(fp)%int = gint
              fpt(fp)%nsl%x = g(gint) * cx(i) + &
                              (1.D0 - g(gint)) * cx(i+delta_i(gint))
              fpt(fp)%nsl%y = g(gint) * cy(j) + &
                              (1.D0 - g(gint)) * cy(j+delta_j(gint))
              fpt(fp)%nsl%z = cz(k)

            END IF

          END DO
          
          k = ord2d(i,j)
          
          fp = fp + 1 
          fpt(fp)%i  = i
          fpt(fp)%j  = j
          fpt(fp)%k  = k
          fpt(fp)%int = 0
          fpt(fp)%nsl%x = cx(i)
          fpt(fp)%nsl%y = cy(j)
          fpt(fp)%nsl%z = topo2d(i,j)

        END DO
      END DO

      RETURN
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      REAL*8 FUNCTION gamma(i_,j_,k_,index)
      USE volcano_topography, ONLY: cx, cy, cz

      IMPLICIT NONE
      ! ... (1 < index < 8) according to the relative position of
      ! ... the neighbours
      INTEGER, INTENT(IN) :: i_, j_, k_, index

      INTEGER :: di_, dj_
      REAL*8  :: num, den

      di_ = delta_i(index)
      dj_ = delta_j(index)

      num = cz(k_) - topo2d(i_,j_)
      den = topo2d(i_+di_, j_+dj_) - topo2d(i_,j_)

      gamma = num / den
      
      RETURN
      END FUNCTION gamma
!----------------------------------------------------------------------
      INTEGER FUNCTION mingam(gam)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: gam(:)
      REAL*8 :: a, ming
      INTEGER :: n

      a = 0.5D0
      ming = 1
      
      DO n = 1, 8
        IF( ABS(a-gam(n)) < ming ) THEN
          ming = ABS(a-gam(n))
          mingam = n
        END IF
      END DO
      
      RETURN
      END FUNCTION mingam
!----------------------------------------------------------------------
!
      END SUBROUTINE forcing3d
!
!----------------------------------------------------------------------
      SUBROUTINE extrafx3d(fx, nfpx, topo2d)
      USE volcano_topography, ONLY: cx, cy, cz, nextx, nexty
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nfpx
      TYPE(forcing_point), DIMENSION(:), INTENT(INOUT) :: fx
      REAL*8, DIMENSION(:,:), INTENT(IN) :: topo2d
      INTEGER :: ijk, i, j, k
      INTEGER :: npx, n
      REAL*8 :: dist1y, dist2y, alpha, beta
      REAL*8 :: quota1, quota2
      LOGICAL :: others
        
      npx = nfpx

      DO ijk = 1, ntot
!
! ... If the vertical velocity is forced but the horizontal
! ... velocity is not, impose the forcing horizontally with
! ... extrapolation ('int'= -10)
!
        others = forcey(ijk) .OR. forcez(ijk)

        IF (others .AND. (.NOT.forcex(ijk))) THEN

          k = ( ijk - 1 ) / ( nx*ny ) + 1
          j = MOD( ijk - 1, nx*ny) / nx + 1
          i = MOD( MOD( ijk - 1, nx*ny ), nx ) + 1

          npx = npx + 1

          fx(npx)%i = i
          fx(npx)%j = j
          fx(npx)%k = k

          IF (cz(k) >= topo2d(i,j)) THEN

            ! ... extrapolation is done if the point lays
            ! ... above the topography ...
            !

            dist1y = cy(j) - ytop(nexty(j)-1)
            dist2y = ytop(nexty(j)) - cy(j)
            beta   = dist1y/(dist1y+dist2y)

            quota1 = beta * ztop2d(nextx(i),nexty(j)) + &
                     (1.D0 - beta) * ztop2d(nextx(i),nexty(j)-1)

            quota2 = beta * ztop2d(nextx(i)-1,nexty(j)) + &
                     (1.D0 - beta) * ztop2d(nextx(i)-1,nexty(j)-1)

            alpha = (cz(k) - quota2) / (quota1 - quota2)

            fx(npx)%nsl%x = alpha * xtop(nextx(i)) + &
                            (1.D0 - alpha) * xtop(nextx(i)-1)

            fx(npx)%nsl%y = cy(j)
            fx(npx)%nsl%z = cz(k)
            fx(npx)%int = 9

          ELSE
            !
            ! ...  do nothing!
            !
            fx(npx)%int = -9
            fx(npx)%nsl%x = 0.D0
            fx(npx)%nsl%y = 0.D0
            fx(npx)%nsl%z = 0.D0
          END IF

        END IF
      END DO
!
      IF (npx/=nfp) CALL error('extrafx','check nfp',1)
!
      RETURN
      END SUBROUTINE extrafx3d
!----------------------------------------------------------------------
      SUBROUTINE extrafy3d(fy, nfpy, topo2d)
      USE volcano_topography, ONLY: cx, cy, cz, nextx, nexty
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nfpy
      TYPE(forcing_point), DIMENSION(:), INTENT(INOUT) :: fy
      REAL*8, DIMENSION(:,:), INTENT(IN) :: topo2d
      INTEGER :: ijk, i, j, k
      INTEGER :: npy, n
      REAL*8 :: dist1x, dist2x, alpha, beta
      REAL*8 :: quota1, quota2
      LOGICAL :: others
        
      npy = nfpy

      DO ijk = 1, ntot
!
! ... If the vertical velocity is forced but the horizontal
! ... velocity is not, impose the forcing horizontally with
! ... extrapolation ('int'= -10)
!
        others = forcex(ijk) .OR. forcez(ijk)

        IF (others .AND. (.NOT.forcey(ijk))) THEN

          k = ( ijk - 1 ) / ( nx*ny ) + 1
          j = MOD( ijk - 1, nx*ny) / nx + 1
          i = MOD( MOD( ijk - 1, nx*ny ), nx ) + 1

          npy = npy + 1

          fy(npy)%i = i
          fy(npy)%j = j
          fy(npy)%k = k
          
          IF (cz(k) >= topo2d(i,j)) THEN
            ! ... extrapolation is done if the point lays
            ! ... above the topography ...
            !
            dist1x = cx(i) - xtop(nextx(i)-1)
            dist2x = xtop(nextx(i)) - cx(i)
            alpha  = dist1x/(dist1x+dist2x)

            quota1 = alpha * ztop2d(nextx(i),nexty(j)) + &
                     (1.D0 - alpha) * ztop2d(nextx(i)-1,nexty(j))

            quota2 = alpha * ztop2d(nextx(i),nexty(j)-1) + &
                     (1.D0 - alpha) * ztop2d(nextx(i)-1,nexty(j)-1)

            beta = (cz(k) - quota2) / (quota1 - quota2)

            fy(npy)%nsl%y = beta * ytop(nexty(j)) + &
                            (1.D0 - beta) * ytop(nexty(j)-1)

            fy(npy)%nsl%x = cx(i)
            fy(npy)%nsl%z = cz(k)
            fy(npy)%int = 10
          ELSE
            ! ...  do nothing!
            !
            fy(npy)%int = -10
            fy(npy)%nsl%x = 0.D0
            fy(npy)%nsl%y = 0.D0
            fy(npy)%nsl%z = 0.D0
          END IF

        END IF
      END DO
!
      IF (npy/=nfp) CALL error('extrafy','check nfp',1)
!
      RETURN
      END SUBROUTINE extrafy3d
!----------------------------------------------------------------------
      SUBROUTINE extrafz3d(fz, nfpz, topo2d)
      USE volcano_topography, ONLY: cx, cy, cz, next
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nfpz
      TYPE(forcing_point), DIMENSION(:), INTENT(INOUT) :: fz
      REAL*8, DIMENSION(:,:), INTENT(IN) :: topo2d
      INTEGER :: ijk, i, j, k
      INTEGER :: npz
      LOGICAL :: others
        
      npz = nfpz

      DO ijk = 1, ntot

        k = ( ijk - 1 ) / ( nx*ny ) + 1
        j = MOD( ijk - 1, nx*ny) / nx + 1
        i = MOD( MOD( ijk - 1, nx*ny ), nx ) + 1

        others = forcex(ijk) .OR. forcey(ijk)
        
        IF (others .AND. (.NOT.forcez(ijk))) THEN
        
          npz = npz + 1
        
          fz(npz)%i = i
          fz(npz)%j = j
          fz(npz)%k = k

          IF (cz(k) >= topo2d(i,j)) THEN
            fz(npz)%int = 11
            fz(npz)%nsl%x = cx(i)
            fz(npz)%nsl%y = cy(j)
            fz(npz)%nsl%z = topo2d(i,j)
          ELSE
            fz(npz)%int = -11
          END IF
        
        END IF
      END DO
!
      IF (npz/=nfp) CALL error('extrafz','check nfp',1)

      RETURN
      END SUBROUTINE extrafz3d
!----------------------------------------------------------------------
      END MODULE immersed_boundaries
!----------------------------------------------------------------------
