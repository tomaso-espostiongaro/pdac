!----------------------------------------------------------------------
      MODULE grid
!
! ... This module contains all the information about the rectilinear
! ... 2D or 3D mesh used for computation, and the cell types 
!
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: dx, dy, dz
      REAL*8, DIMENSION(:), ALLOCATABLE :: indx, indy, indz
      REAL*8, DIMENSION(:), ALLOCATABLE :: r, rb
      REAL*8, DIMENSION(:), ALLOCATABLE :: inr, inrb
      REAL*8, DIMENSION(:), ALLOCATABLE :: x, xb
      REAL*8, DIMENSION(:), ALLOCATABLE :: y, yb
      REAL*8, DIMENSION(:), ALLOCATABLE :: z, zb
!
      INTEGER :: itc

      TYPE blbody
        INTEGER :: xlo
        INTEGER :: xhi
        INTEGER :: ylo
        INTEGER :: yhi
        INTEGER :: zlo
        INTEGER :: zhi
        INTEGER :: typ
      END TYPE
!
      TYPE (blbody), ALLOCATABLE :: iob(:)
!
      INTEGER, DIMENSION(:), ALLOCATABLE :: fl   ! temporary global array
      INTEGER, DIMENSION(:), ALLOCATABLE :: flag ! local flag array
!
      INTEGER :: west, east, south, north, bottom, top
      CHARACTER(LEN=80) :: topography

      REAL*8 :: zzero
!
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_grid
!
      USE dimensions, ONLY: nx, ny, nz, ntot
      USE control_flags, ONLY: job_type
!
      IMPLICIT NONE
!
! ...   Set the appropriate total number of cells
! ...   and the coordinate system

      IF( job_type == '2D' ) THEN
        ntot = nx*nz
      ELSE IF( job_type == '3D' ) THEN
        ntot = nx*ny*nz
      ELSE
        CALL error( ' allocate_grid ', ' wrong job_type '//job_type, 1)
      END IF

      ALLOCATE( dx(nx), dy(ny), dz(nz) ) 
      ALLOCATE( indx(nx), indy(ny), indz(nz) )
      ALLOCATE( r(nx), rb(nx) )
      ALLOCATE( inr(nx), inrb(nx) )
      ALLOCATE( x(nx), xb(nx) )
      ALLOCATE( y(ny), yb(ny) )
      ALLOCATE( z(nz), zb(nz) )
      ALLOCATE( fl(ntot) )

      RETURN
      END SUBROUTINE allocate_grid
!----------------------------------------------------------------------
      SUBROUTINE allocate_blbody
      USE dimensions
      IMPLICIT NONE
      
      ALLOCATE(iob(no))
      
      RETURN
      END SUBROUTINE allocate_blbody
!----------------------------------------------------------
      SUBROUTINE grid_setup
!
      USE control_flags, ONLY: job_type
      USE dimensions, ONLY: nx, ny, nz

      REAL*8, PARAMETER :: VERYBIG = 1.0d+10
      INTEGER :: i, j, k, ijk
      REAL*8 :: zrif

      OPEN(17,FILE='mesh.dat')
!
! ... Compute grid point (x,y,z) locations
!
      xb(1)=0.D0
      x(1)=-0.5D0*dx(1)+xb(1)
      DO i=2,nx
        xb(i)=(xb(i-1)+dx(i))
        x(i)=(xb(i)-0.5D0*dx(i))
      END DO

      yb(1) = 0.D0
      y(1)=-0.5D0*dy(1)+yb(1)
      DO j=2,ny
        yb(j) = yb(j-1) + dy(j)
        y(j)=(yb(j)-0.5D0*dy(j))
      END DO

      zb(1) = zzero
      z(1)=-0.5D0*dz(1)+zb(1)
      DO k=2,nz
        zb(k) = zb(k-1) + dz(k)
        z(k)=(zb(k)-0.5D0*dz(k))
      END DO

      WRITE(17,*) x
      WRITE(17,*) y
      WRITE(17,*) z
      WRITE(17,*) xb
      WRITE(17,*) yb
      WRITE(17,*) zb

      CLOSE(17)

! ... inverse of the cell sizes

      indx = 1.D0 / dx
      indz = 1.D0 / dz

      IF (job_type == '3D') THEN
        indy = 1.D0 / dy
      ELSE IF (job_type == '2D') THEN
        indy = 0.D0
      END IF
!
! ... Set "r" and "rb" absolute coordinate
! ... (only for 2D cylindrical coordinates)
!
      IF((itc == 0) .OR. (job_type == '3D')) THEN

        r(:) = 1.D0
        rb(:) = 1.D0
        inr(:) = 1.D0
        inrb(:) = 1.D0

      ELSE IF ((itc == 1) .AND. (job_type == '2D')) THEN

        r(:)  = x(:)
        rb(:) = xb(:)

        IF (rb(1) == 0.0D0) THEN
          inrb(1) = VERYBIG
        ELSE
          inrb(1) = 1.0D0/rb(1)
        END IF
        DO i=2,nx
          inrb(i)=1.D0/rb(i)
        END DO  
        DO i=1,nx
          inr(i)=1.D0/r(i)
        END DO  

      END IF
!
! ... Set cell flags
!
      CALL flic
!      
      RETURN 
      END SUBROUTINE grid_setup
!----------------------------------------------------------------------
      SUBROUTINE flic
      USE dimensions, ONLY: nx, ny, nz, no
      USE control_flags, ONLY: job_type
!
      IMPLICIT NONE
      INTEGER :: i, j, k, n, ijk
!
! ... Default cell type corresponds to fluid flow
!
      fl = 1
!
      IF( job_type == '2D' ) THEN
        !
        ! ... Specify cell flags on mesh boundaries
        !
        DO i = 1, nx
          k = 1
          ijk = i + (k-1) * nx
          fl(ijk) = bottom
          k = nz
          ijk = i + (k-1) * nx
          fl(ijk) = top
        END DO
        !
        DO k = 1, nz
          i = 1
          ijk = i + (k-1) * nx
          fl(ijk) = west
          i = nx
          ijk = i + (k-1) * nx
          fl(ijk) = east
        END DO
        !
        ! ... Specify cell flags for blocks
        !
        DO n = 1, no
          DO k = iob(n)%zlo, iob(n)%zhi 
            DO i = iob(n)%xlo, iob(n)%xhi
              ijk = i + (k-1) * nx

              fl(ijk) = iob(n)%typ
              
            END DO
          END DO
        END DO
      ELSE IF( job_type == '3D' ) THEN
        !
        ! ... Specify cell flags on mesh boundaries
        !
        DO i = 1, nx
          DO j = 1, ny
            k = 1
            ijk = i + ( (j-1) + (k-1)*ny ) * nx
            fl(ijk) = bottom
            k = nz
            ijk = i + ( (j-1) + (k-1)*ny ) * nx
            fl(ijk) = top
          END DO
        END DO
        !
        DO j = 1, ny
          DO k = 1, nz
            i = 1
            ijk = i + ( (j-1) + (k-1)*ny ) * nx
            fl(ijk) = west
            i = nx
            ijk = i + ( (j-1) + (k-1)*ny ) * nx
            fl(ijk) = east
          END DO
        END DO
        !
        DO k = 1, nz
          DO i = 1, nx
            j = 1
            ijk = i + ( (j-1) + (k-1)*ny ) * nx
            fl(ijk) = south
            j = ny
            ijk = i + ( (j-1) + (k-1)*ny ) * nx
            fl(ijk) = north
          END DO
        END DO
        !
        ! ... Specify cell flags for blocks
        !
        DO n = 1, no
          DO k = iob(n)%zlo, iob(n)%zhi
            DO j = iob(n)%ylo, iob(n)%yhi
              DO i = iob(n)%xlo, iob(n)%xhi
                ijk = i + ( (j-1) + (k-1)*ny ) * nx

                fl(ijk) = iob(n)%typ
              
              END DO
            END DO
          END DO
        END DO
        !
      END IF
!
! ... Set flags on corners
!
      IF( job_type == '2D' ) THEN
        IF ( fl(nz*nx - 1) == 4 .AND. fl((nz-1)*nx) == 4) THEN
          fl(nz*nx) = 4
        END IF
      END IF
      
      RETURN
      END SUBROUTINE flic
!----------------------------------------------------------------------
      END MODULE grid
!----------------------------------------------------------------------
