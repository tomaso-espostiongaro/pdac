!----------------------------------------------------------------------
      MODULE immersed_boundaries
!----------------------------------------------------------------------
      USE dimensions, ONLY: nx, ny, nz, ntot, ntr
      USE grid, ONLY: fl
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
! ... In 2D: -3=lin sx/top, -2=linsx, -1=bilsx, 
! ...        0=lintop, 1=bildx, 2=lindx, 3=lin dx/top 
! ... In 3D: always linear interpolation with the nerighbour
! ...        (from 1 to 8) that is more convenient       
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
!
! ... This array associates to each local mesh point an integer that
! ... is the decimal representation of the binary number giving the
! ... number of cell faces laying outside the topography
      INTEGER*2, ALLOCATABLE :: bd(:)
      REAL*8, ALLOCATABLE :: bdr(:,:)
      REAL*8, ALLOCATABLE :: vf(:)

      INTEGER :: immb           
!
      INTERFACE faces
        MODULE PROCEDURE faces_real, faces_int
      END INTERFACE
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE set_forcing
      USE control_flags, ONLY: job_type, lpr
      USE volcano_topography, ONLY: interpolate_2d, interpolate_dem
      USE volcano_topography, ONLY: ord, next
      USE volcano_topography, ONLY: ord2d, nextx, nexty
      USE grid, ONLY: x, xb, y, yb, z, zb

      IMPLICIT NONE
      INTEGER :: p, np
      INTEGER :: i,j,k,ijk
      INTEGER :: nfpx, nfpy, nfpz
!
! ... Conditional array to be used in loops:
! ... this value is .TRUE. when force has to be computed at a given location
!
      LOGICAL, ALLOCATABLE :: force(:)
!
! ... If Immersed Boundaries are used, identify the forcing points
! ... and set interpolation parameters
!
      !
      ! ... Allocate the logical array that is used to 
      ! ... identify the forcing points
      !
      ALLOCATE(force(ntot)); force = .FALSE.

      IF (job_type == '2D') THEN

        ALLOCATE(topo_c(nx))
        ALLOCATE(topo_x(nx))

        !
        ! ... interpolate the topography on x-staggered mesh
        ! ... and count the forcing points along x
        !
        CALL interpolate_2d(xb, z, topo_x, force)
        nfpx = COUNT(force)
        ALLOCATE(fptx(nfpx))
        !
        ! ... Forcing along x
        CALL forcing2d(xb, z, topo_x, fptx)
        !
        ! ... Set flag = 1 on forcing points
        DO np = 1, nfpx
          i = fptx(np)%i
          k = fptx(np)%k
          ijk = i + (k-1) * nx
          IF (k>1) fl(ijk) = 1
        END DO

        !
        ! ... interpolate the topography on z-staggered mesh
        ! ... and count the forcing points along z
        !
        CALL interpolate_2d(x, zb, topo_c, force)
        nfpz = COUNT(force)
        ALLOCATE(fptz(nfpz))
        !
        ! ... Forcing along z
        CALL forcing2d(x, zb, topo_c, fptz)
        !
        ! ... Set flag = 1 on forcing points
        DO np = 1, nfpz
          i = fptz(np)%i
          k = fptz(np)%k
          ijk = i + (k-1) * nx
          IF (k>1) fl(ijk) = 1
        END DO
        !
      ELSE IF (job_type == '3D') THEN

        ALLOCATE(topo2d_c(nx,ny))
        ALLOCATE(topo2d_x(nx,ny))
        ALLOCATE(topo2d_y(nx,ny))
        !
        ! ... interpolate the topography on x-staggered mesh
        ! ... and count the forcing points along x
        !
        CALL interpolate_dem(xb, y, z, topo2d_x, force)
        nfpx = COUNT(force)
        ALLOCATE(fptx(nfpx))
        !
        ! ... Forcing along x
        CALL forcing3d(xb, y, z, topo2d_x, fptx)
        !
        ! ... Set flag = 1 on forcing points
        DO np = 1, nfpx
          i = fptx(np)%i
          j = fptx(np)%j
          k = fptx(np)%k
          ijk = i + (j-1) * nx + (k-1) * nx * ny
          IF (k>1) fl(ijk) = 1
        END DO
        
        !
        ! ... interpolate the topography on y-staggered mesh
        ! ... and count the forcing points along y
        !
        CALL interpolate_dem(x, yb, z, topo2d_y, force)
        nfpy = COUNT(force)
        ALLOCATE(fpty(nfpy))
        !
        ! ... Forcing along y
        CALL forcing3d(x, yb, z, topo2d_y, fpty)
        !
        ! ... Set flag = 1 on forcing points
        DO np = 1, nfpy
          i = fpty(np)%i
          j = fpty(np)%j
          k = fpty(np)%k
          ijk = i + (j-1) * nx + (k-1) * nx * ny
          IF (k>1) fl(ijk) = 1
        END DO
        
        !
        ! ... interpolate the topography on z-staggered mesh
        ! ... and count the forcing points along z
        !
        CALL interpolate_dem(x, y, zb, topo2d_c, force)
        nfpz = COUNT(force)
        ALLOCATE(fptz(nfpz))
        !
        ! ... Forcing along z
        CALL forcing3d(x, y, zb, topo2d_c, fptz)
        !
        ! ... Set flag = 1 on forcing points
        DO np = 1, nfpz
          i = fptz(np)%i
          j = fptz(np)%j
          k = fptz(np)%k
          ijk = i + (j-1) * nx + (k-1) * nx * ny
          IF (k>1) fl(ijk) = 1
        END DO
        
      END IF

      IF (lpr > 1) THEN
        IF (mpime == root) THEN
          OPEN(UNIT=15,FILE='fptx.dat',STATUS='UNKNOWN')
          IF (job_type == '3D') &
            OPEN(UNIT=16,FILE='fpty.dat',STATUS='UNKNOWN')
          OPEN(UNIT=17,FILE='fptz.dat',STATUS='UNKNOWN')
          DO np = 1, SIZE(fptx)
            WRITE(15,33) np, fptx(np)
          END DO
          CLOSE(15)
          IF (job_type == '3D') THEN
            DO np = 1, SIZE(fpty)
              WRITE(16,33) np, fpty(np)
            END DO
            CLOSE(16)
          END IF
          DO np = 1, SIZE(fptz)
            WRITE(17,33) np, fptz(np)
          END DO
          CLOSE(17)
        END IF
 33   FORMAT(5(I6),4(F18.3))
      END IF
!
      DEALLOCATE (force)
!
      RETURN
      END SUBROUTINE set_forcing
!----------------------------------------------------------------------
      SUBROUTINE forcing2d(cx, cz, topo, fpt)
      USE volcano_topography, ONLY: ord, next
!
! ... locate the forcing points and the type of interpolation
! ... locate the points where no-slip condition is applied

      IMPLICIT NONE

      REAL*8, DIMENSION(:), INTENT(IN) :: cx, cz
      REAL*8, DIMENSION(:), INTENT(IN) :: topo
      TYPE(forcing_point), DIMENSION(:), INTENT(OUT) :: fpt

      INTEGER :: i, fp
      REAL*8 :: t,s
!
! ... Skip the -- first grid point -- 
! ... Loop over grid points in x-direction
!
      fp = 0
      DO i=2,nx-1
         
         IF ((ord(i)-ord(i-1) > 0).AND.(topo(i+1)-topo(i-1) >= 0)) THEN
            
            CALL increasing_profile(i,fp,fpt)
            
         ELSEIF ( ord(i+1)-ord(i) < 0) THEN

               CALL decreasing_profile(i,fp,fpt)

         ELSE 
            
            CALL stationnary_profile(i,fp,fpt)
            
         END IF
         
      END DO
!     
! ... Skip the -- last grid point -- 
!
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
        IF (err > 0) WRITE(7,*) 'Error in fp: ', fp

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
        IF (err > 0) WRITE(7,*) 'Error in fp: ', fp

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
        IF (err > 0) WRITE(7,*) 'Error in fp: ', fp

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
        IF (err > 0) WRITE(7,*) 'Error in fp: ', fp

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
!
      END SUBROUTINE forcing2d
!
!----------------------------------------------------------------------
      SUBROUTINE forcing3d(cx, cy, cz, topo2d, fpt)
      USE volcano_topography, ONLY: ord2d, next
!
! ... locate the forcing points and the type of interpolation
! ... locate the points where no-slip condition is applied
!
      IMPLICIT NONE

      REAL*8, DIMENSION(:), INTENT(IN) :: cx, cy, cz
      REAL*8, DIMENSION(:,:), INTENT(IN) :: topo2d
      TYPE(forcing_point), DIMENSION(:), INTENT(OUT) :: fpt

      INTEGER :: i, j, k, fp
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

! ... Skip the boundary points
!
      DO i = 2, nx - 1
        DO j = 2, ny - 1
          DO k = 2, ord2d(i,j) - 1
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
      SUBROUTINE faces_int(ijk, b_e, b_w, b_t, b_b, b_n, b_s, ivf)
      USE control_flags, ONLY: job_type
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: ijk
      INTEGER, INTENT(OUT) :: b_e, b_w, b_t, b_b, b_n, b_s
      REAL*8, INTENT(OUT) :: ivf
      INTEGER*2 :: num
 
      b_e = 0
      num = 1
      IF( IAND(bd(ijk),num) /= 0 )  b_e = 1 
 
      b_w = 0
      num = 2
      IF( IAND(bd(ijk),num) /= 0 )  b_w = 1
 
      b_t = 0
      num = 4
      IF( IAND(bd(ijk),num) /= 0 )  b_t = 1
 
      b_b = 0
      num = 8
      IF( IAND(bd(ijk),num) /= 0 )  b_b = 1
      
      b_n = 0
      num = 16
      IF( IAND(bd(ijk),num) /= 0 ) b_n = 1
 
      b_s = 0
      num = 32
      IF( IAND(bd(ijk),num) /= 0 ) b_s = 1
 
      IF (vf(ijk) > 0.D0) THEN
        ivf = 1.D0 / vf(ijk)
      ELSE
        ivf = 0.D0
      END IF
 
      RETURN
      END SUBROUTINE faces_int
!----------------------------------------------------------------------
      SUBROUTINE faces_real(ijk, b_e, b_w, b_t, b_b, b_n, b_s, ivf)
      USE control_flags, ONLY: job_type
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: ijk
      REAL*8, INTENT(OUT) :: b_e, b_w, b_t, b_b, b_n, b_s
      REAL*8, INTENT(OUT) :: ivf
 
      b_e = bdr(ijk,1)
      b_w = bdr(ijk,2)
      b_t = bdr(ijk,3)
      b_b = bdr(ijk,4)
      b_n = bdr(ijk,5)
      b_s = bdr(ijk,6)

      IF (vf(ijk) > 0.D0) THEN
        ivf = 1.D0 / vf(ijk)
      ELSE
        ivf = 0.D0
      END IF
 
      RETURN
      END SUBROUTINE faces_real
!----------------------------------------------------------------------
      SUBROUTINE faces_real_2(ijk, b_e, b_w, b_t, b_b, b_n, b_s, ivf)
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: ijk
      REAL*8, INTENT(OUT) :: b_e, b_w, b_t, b_b, b_n, b_s
      REAL*8, INTENT(OUT) :: ivf
 
      ivf = 0.D0
 
      b_e = bdr(ijk,1)
      IF (b_e  <= 0.D0) THEN
        b_e = b_e + DSQRT( b_e**2 - b_e )
      ELSE
        b_e = b_e - DSQRT( b_e**2 - b_e )
      END IF
      b_e = b_e / ( 2.D0 * b_e - 1.D0 )

      b_w = bdr(ijk,2)
      IF (b_w  <= 0.D0) THEN
        b_w = b_w + DSQRT( b_w**2 - b_w )
      ELSE
        b_w = b_w - DSQRT( b_w**2 - b_w )
      END IF
      b_w = b_w / ( 2.D0 * b_w - 1.D0 )

      b_t = bdr(ijk,3)
      IF (b_t  <= 0.D0) THEN
        b_t = b_t + DSQRT( b_t**2 - b_t )
      ELSE
        b_t = b_t - DSQRT( b_t**2 - b_t )
      END IF
      b_t = b_t / ( 2.D0 * b_t - 1.D0 )

      b_b = bdr(ijk,4)
      IF (b_b  <= 0.D0) THEN
        b_b = b_b + DSQRT( b_b**2 - b_b )
      ELSE
        b_b = b_b - DSQRT( b_b**2 - b_b )
      END IF
      b_b = b_b / ( 2.D0 * b_b - 1.D0 )

      b_n = bdr(ijk,5)
      IF (b_n  <= 0.D0) THEN
        b_n = b_n + DSQRT( b_n**2 - b_n )
      ELSE
        b_n = b_n - DSQRT( b_n**2 - b_n )
      END IF
      b_n = b_n / ( 2.D0 * b_n - 1.D0 )

      b_s = bdr(ijk,2)
      IF (b_s  <= 0.D0) THEN
        b_s = b_s + DSQRT( b_s**2 - b_s )
      ELSE
        b_s = b_s - DSQRT( b_s**2 - b_s )
      END IF
      b_s = b_s / ( 2.D0 * b_s - 1.D0 )

      RETURN
      END SUBROUTINE faces_real_2
!----------------------------------------------------------------------
      END MODULE immersed_boundaries
!----------------------------------------------------------------------
