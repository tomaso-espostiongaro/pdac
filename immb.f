!----------------------------------------------------------------------
      MODULE immersed_boundaries
!----------------------------------------------------------------------
      USE dimensions, ONLY: nx, ny, nz, ntot, ntr
      USE parallel, ONLY: mpime, root
      USE volcano_topography, ONLY: xtop, ytop, ztop
      USE volcano_topography, ONLY: cx, cy, cz

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
! ... to be done (-10= extrap z -3=lin sx/top, -2=linsx, -1=bilsx, 
! ... 0=lintop, 1=bildx, 2=lindx, 3=lin dx/top, 10=extrap x), 
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

      TYPE(forcing_point), ALLOCATABLE :: fptx(:), fptz(:)
      TYPE(forcing_point), ALLOCATABLE :: forx(:), forz(:)
!
! ... Conditional array to be used in loops:
! ... this value is .TRUE. when force has to be computed at a given location
!
      LOGICAL, ALLOCATABLE :: forcex(:)
      LOGICAL, ALLOCATABLE :: forcez(:)
      LOGICAL, ALLOCATABLE :: forced(:)
!
! ... This function gives the increasing number of local forcing points 
      INTEGER, ALLOCATABLE :: numx(:), numz(:)
      INTEGER :: nfp
!
      REAL*8, ALLOCATABLE  :: topo_c(:)  ! topography on the grid
      REAL*8, ALLOCATABLE  :: topo_x(:)  ! topography on the grid
!           
      REAL*8, ALLOCATABLE  :: topo2d_c(:,:)  ! topography on the grid
      REAL*8, ALLOCATABLE  :: topo2d_x(:,:)  ! topography on the grid
      REAL*8, ALLOCATABLE  :: topo2d_y(:,:)  ! topography on the grid
!           
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE import_topo
      USE control_flags, ONLY: job_type, lpr, immb
      USE volcano_topography, ONLY: grid_locations, dist
      USE volcano_topography, ONLY: interpolate_2d, ord, next
      USE volcano_topography, ONLY: interpolate_dem, ord2d, nextx, nexty

      IMPLICIT NONE
      INTEGER :: p, nfpx, nfpz
      INTEGER :: i,j,k,ijk
!
      ALLOCATE(dist(ntot))
!
! ... interpolate topography on cell centers and write the
! ... implicit profile
!
      IF (job_type == '2D') THEN
  
        ALLOCATE(topo_c(nx))
        ALLOCATE(topo_x(nx))

        ALLOCATE(next(nx))
        ALLOCATE(ord(nx))

        CALL grid_locations(0,0,0)
        CALL interpolate_2d(topo_c)

      ELSE IF (job_type == '3D') THEN

        ALLOCATE(topo2d_c(nx,ny))
        ALLOCATE(topo2d_x(nx,ny))
        ALLOCATE(topo2d_x(nx,ny))

        ALLOCATE(nextx(nx))
        ALLOCATE(nexty(ny))
        ALLOCATE(ord2d(nx,ny))

        CALL grid_locations(0,0,0)
        CALL interpolate_dem(topo2d_c)

        WRITE(6,*) 'TOPOGRAPHY'
        DO j= 1, ny
          DO i= 1, nx
            WRITE(6,*) i, j, topo2d_c(i,j)
          END DO
        END DO

      END IF

      IF (mpime == root) THEN
        OPEN(UNIT=14,FILE='improfile.dat',STATUS='UNKNOWN')
        WRITE(14,*) dist
        CLOSE(14)
      END IF
!
      CALL error('immb','debug',1)
      IF (immb >= 1) THEN
!
! ... These arrays are exported and scattered among processors
! ... after domain decomposition
!
        ALLOCATE(forcex(ntot))
        ALLOCATE(forcez(ntot))
!
! ... Interpolate topography on staggered points
! ... If immersed boundary technique is used, locate forcing points
! ... and set interpolation parameters
!
        !
        ! ... interpolate the topography on x-staggered mesh
        CALL grid_locations(1,0,0)
        CALL interpolate_2d(topo_x,forcex)
        nfpx = COUNT(forcex)
        ALLOCATE(fptx(nfpx))
        !
        ! ... Forcing along x
        CALL forcing(fptx,topo_x)
! 
        !
        ! ... interpolate the topography on z-staggered mesh
        CALL grid_locations(0,0,1)
        CALL interpolate_2d(topo_c,forcez)
        nfpz = COUNT(forcez)
        ALLOCATE(fptz(nfpz))
        !
        ! ... Forcing along z
        CALL forcing(fptz,topo_c)
!
! ... Merge forcing points subsets
!
        nfp = COUNT( forcex .OR. forcez )
        ALLOCATE (forx (nfp) )
        ALLOCATE (forz (nfp) )

        CALL grid_locations(1,0,0)
        CALL interpolate_2d(topo_x)
        forx(1:nfpx) = fptx(1:nfpx)
        CALL extrafx(forx, nfpx)

        CALL grid_locations(0,0,1)
        CALL interpolate_2d(topo_c)
        forz(1:nfpz) = fptz(1:nfpz)
        CALL extrafz(forz, nfpz)

        IF (mpime == root) THEN
          OPEN(UNIT=15,FILE='forx.dat',STATUS='UNKNOWN')
          OPEN(UNIT=16,FILE='forz.dat',STATUS='UNKNOWN')
          DO p = 1, nfp
            WRITE(15,33) p, forx(p)
            WRITE(16,33) p, forz(p)
          END DO
          CLOSE(15)
          CLOSE(16)
 33       FORMAT(5(I5),4(F13.6))
        END IF

        DEALLOCATE(fptx)
        DEALLOCATE(fptz)

      END IF
! 
! ... Set flags to identify the topography (where transport equations
! ... are not solved )
!
      CALL set_flag3(topo_c)
!
      DEALLOCATE (xtop, ztop)
      DEALLOCATE (cx, cy, cz)
      DEALLOCATE (next, ord, dist)
!
      RETURN
      END SUBROUTINE import_topo
!----------------------------------------------------------------------
      SUBROUTINE forcing(fpt, topo)
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
      END SUBROUTINE forcing
!----------------------------------------------------------------------
      SUBROUTINE extrafx(fx, nfpx)
      USE volcano_topography, ONLY: cx, cy, cz, next
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nfpx
      TYPE(forcing_point), DIMENSION(:), INTENT(INOUT) :: fx
      INTEGER :: ijk, i, j, k
      INTEGER :: npx, n
      REAL*8 :: grad
        
      npx = nfpx

      DO ijk = 1, ntot
        k = ( ijk - 1 ) / nx + 1
        i = MOD( ( ijk - 1 ), nx) + 1
        IF (forcez(ijk).AND.(.NOT.forcex(ijk))) THEN
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

          IF (cz(k) >= topo_x(i)) THEN
            fx(npx)%int = -10
          ELSE
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
      END SUBROUTINE extrafx
!----------------------------------------------------------------------
      SUBROUTINE extrafz(fz, nfpz)
      USE volcano_topography, ONLY: cx, cy, cz, next
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nfpz
      TYPE(forcing_point), DIMENSION(:), INTENT(INOUT) :: fz
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
          IF (cz(k) >= topo_c(i)) THEN
            fz(npz)%int = 10
            fz(npz)%nsl%x = cx(i)
            fz(npz)%nsl%z = topo_c(i)
          ELSE
            fz(npz)%int = 100
          END IF
        END IF
      END DO
!
      IF (npz/=nfp) CALL error('extrafz','check nfp',1)

      RETURN
      END SUBROUTINE extrafz
!----------------------------------------------------------------------
      SUBROUTINE set_flag3(topo)
!
! ... Set cell-flag = 3 in those cells belonging to the topography
! ... where forcing is not applied. Values of fields in these cells
! ... are set to zero when initialized or kept undefined
!
      USE control_flags, ONLY: lpr, immb
      USE grid, ONLY: fl, z
      IMPLICIT NONE

      REAL*8, DIMENSION(:), INTENT(IN) :: topo
      INTEGER :: i, k, ijk
      LOGICAL :: undg
!
! ... If immersed boundaries are used, set flag = 3 below forcing points
! ... otherwise set flag = 3 when topography is above the cell center
! ... (Assumes topography on domain bottom!)
!
      IF (immb >= 1) THEN
        DO i=2, nx-1
          DO k = 1, nz
            ijk = i + (k-1) * nx
            IF ((.NOT.forcex(ijk)).AND.(.NOT.forcez(ijk))) THEN
              fl(ijk) = 3
            ELSE
              EXIT
            END IF
          END DO
        END DO
      ELSE IF (immb < 1) THEN
        DO i=2, nx-1
          DO k = 1, nz
            ijk = i + (k-1) * nx
            IF (topo(i) > z(k)) THEN
              fl(ijk) = 3
            ELSE
              EXIT
            END IF
          END DO
        END DO
      END IF
!
      IF (lpr > 1 .AND. (mpime == root)) THEN
        DO k = nz, 1, -1
          WRITE(6,11) (fl(i + (k-1) * nx), i=1,nx)
        END DO
      END IF
                                                                               
 11   FORMAT(80(I1))

      END SUBROUTINE set_flag3
!----------------------------------------------------------------------
      END MODULE immersed_boundaries
!----------------------------------------------------------------------
