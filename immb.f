!----------------------------------------------------------------------
      MODULE immersed_boundaries
!----------------------------------------------------------------------
      USE dimensions, ONLY: nx, ny, nz, ntot
      USE parallel, ONLY: mpime, root

      IMPLICIT NONE

      TYPE point
        REAL*8 :: x
        REAL*8 :: y
        REAL*8 :: z
      END TYPE point
!
! ... Assemble parameters identifying a forcing point: (i,j,k) are the 
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
! ... Parameters identifying the topography interpolated on the
! ... numerical mesh.
!
      REAL*8, ALLOCATABLE, DIMENSION(:) :: cx, cy, cz
      
      REAL*8, ALLOCATABLE  :: topo_c(:)  ! topography on the grid
      REAL*8, ALLOCATABLE  :: topo_x(:)  ! topography on the grid
      INTEGER, ALLOCATABLE :: next(:)  ! index of the topographic point
                                       ! located right of a grid point 
      INTEGER, ALLOCATABLE :: ord(:)   ! z-coord of grid locations laying
                                       ! below the profile
      REAL*8, ALLOCATABLE  :: dist(:)  ! vertical distance above topography

      REAL*8, ALLOCATABLE :: xtop(:),ztop(:) ! input topography points
      INTEGER :: noditop
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
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE import_topo
      USE control_flags, ONLY: job_type, lpr, immb

      IMPLICIT NONE
      INTEGER :: p, nfpx, nfpz
      INTEGER :: i,j,k,ijk
      LOGICAL, ALLOCATABLE :: temp(:)
!
      IF (job_type == '3D') RETURN
!
! ... Read the topography file (at arbitrary resolution)
!
      CALL read_profile

      ALLOCATE(next(nx))
      ALLOCATE(ord(nx))
      ALLOCATE(topo_c(nx))
      ALLOCATE(topo_x(nx))
      ALLOCATE(dist(ntot))
      ALLOCATE(temp(ntot))
!
! ... interpolate topography on cell centers and write the
! ... implicit profile
!
      CALL grid_locations(0,0,0)
      CALL intertop(topo_c,temp)
      
      IF (mpime == root) THEN
        OPEN(UNIT=14,FILE='improfile.dat',STATUS='UNKNOWN')
        WRITE(14,*) dist
        CLOSE(14)
      END IF

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
! ... Forcing along x
!
        CALL grid_locations(1,0,0)
        CALL intertop(topo_x,forcex)
        nfpx = COUNT(forcex)
        ALLOCATE(fptx(nfpx))
        CALL forcing(fptx,topo_x)
! 
! ... Forcing along z
!
        CALL grid_locations(0,0,1)
        CALL intertop(topo_c,forcez)
        nfpz = COUNT(forcez)
        ALLOCATE(fptz(nfpz))
        CALL forcing(fptz,topo_c)
!
! ... Merge forcing points subsets
!
        nfp = COUNT( forcex .OR. forcez )
        ALLOCATE (forx (nfp) )
        ALLOCATE (forz (nfp) )
        forx(1:nfpx) = fptx(1:nfpx)
        forz(1:nfpz) = fptz(1:nfpz)
        CALL extraf(forx, forz, nfpx, nfpz, topo_x, topo_c)

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

      DEALLOCATE (xtop, ztop)
      DEALLOCATE (cx, cy, cz)
      DEALLOCATE (next, ord, dist, temp)
! 
! ... Set flags to identify the topography (where transport equations
! ... are not solved )
!
      CALL set_flag3(topo_c)

      RETURN
      END SUBROUTINE import_topo
!----------------------------------------------------------------------
      SUBROUTINE read_profile
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

      END SUBROUTINE read_profile
!----------------------------------------------------------------------
      SUBROUTINE intertop(topo,ff)
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
      END SUBROUTINE intertop
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
! ... Set flag 3 below specified flow cells
!
      DO i=1, nx
        undg = .FALSE.
        DO k = nz, 1, -1
          ijk = i + (k-1) * nx
          IF (fl(ijk)/=5 .AND. undg) fl(ijk) = 3
          IF ( fl(ijk) == 5 ) undg = .TRUE.
        END DO
      END DO
!
      IF (lpr > 1 .AND. (mpime == root)) THEN
        DO k = nz, 1, -1
          WRITE(6,11) (fl(i + (k-1) * nx), i=1,nx)
        END DO
      END IF
                                                                               
 11   FORMAT(80(I1))

      END SUBROUTINE set_flag3
!----------------------------------------------------------------------
      SUBROUTINE forcing(fpt, topo)
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
      SUBROUTINE extraf(fx, fz, nfpx, nfpz, topx, topc)
!
      USE grid, ONLY: x, xb, z, zb, dx, dz

      IMPLICIT NONE
      REAL*8, DIMENSION(:), INTENT(IN) :: topx, topc
      INTEGER, INTENT(IN) :: nfpx, nfpz
      TYPE(forcing_point), DIMENSION(:), INTENT(INOUT) :: fx, fz
      INTEGER :: ijk, i, j, k
      INTEGER :: npx, npz
      REAL*8 :: gradtop, deltop
        
      npx = nfpx
      npz = nfpz

      DO ijk = 1, ntot
        ! k = ( ijk - 1 ) / ( nx*ny ) + 1
        ! j = MOD( ijk - 1, nx*ny) / nx + 1
        ! i = MOD( MOD( ijk - 1, nx*ny ), nx ) + 1
        k = ( ijk - 1 ) / nx + 1
        i = MOD( ( ijk - 1 ), nx) + 1
        IF (forcex(ijk).AND.(.NOT.forcez(ijk))) THEN
          npz = npz + 1
          fz(npz)%i = i
          fz(npz)%k = k
          IF (zb(k) >= topc(i)) THEN
            fz(npz)%int = 10
            fz(npz)%nsl%x = x(i)
            fz(npz)%nsl%z = topc(i)
          ELSE
            fz(npz)%int = 100
          END IF
        ELSE IF (forcez(ijk).AND.(.NOT.forcex(ijk))) THEN
          npx = npx + 1
          fx(npx)%i = i
          fx(npx)%k = k
          gradtop = ( topx(i) - topc(i) ) / ( 0.5D0 * dx(i) ) 
          deltop  = ( topx(i) - z(k) )
          IF (z(k) >= topx(i)) THEN
            fx(npx)%int = -10
            fx(npx)%nsl%x = xb(i) - deltop / gradtop
            fx(npx)%nsl%z = z(k)
          ELSE
            fx(npx)%int = -100
          END IF
        END IF
      END DO
      IF (npx/=nfp.OR.npz/=nfp) CALL error('extraf','check nfp',1)

      RETURN
      END SUBROUTINE extraf
!----------------------------------------------------------------------
      SUBROUTINE grid_locations(sx,sy,sz)
      USE dimensions, ONLY: nx, ny, nz
      USE grid, ONLY: x, y, z, xb, yb, zb

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
      END MODULE immersed_boundaries
!----------------------------------------------------------------------
