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
        REAL*8, DIMENSION(:), ALLOCATABLE :: x, xb, zb
        REAL*8, DIMENSION(:), ALLOCATABLE :: inx, inxb
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
!*******MX3D fl da eliminare in un secondo tempo
!
        INTEGER, DIMENSION(:), ALLOCATABLE :: fl

        INTEGER, DIMENSION(:), ALLOCATABLE :: fl_l
!
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------

      SUBROUTINE allocate_grid
!
        USE dimensions
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
        ALLOCATE( x(nx), xb(nx) )
        ALLOCATE( inx(nx), inxb(nx) )
        ALLOCATE( zb(nz) )
        ALLOCATE( fl(ntot) )

        RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
!
      SUBROUTINE allocate_blbody
        USE dimensions
        IMPLICIT NONE
          ALLOCATE(iob(no))
        RETURN
      END SUBROUTINE
!
!----------------------------------------------------------
!
        SUBROUTINE grid_setup(zzero)
!
          USE dimensions
          USE control_flags, ONLY: job_type

          REAL*8, INTENT(IN) :: zzero
          REAL*8, PARAMETER :: VERYBIG = 1.0d+10
          INTEGER :: i, j, k
!
! ... Set "x" and "xb" absolute coordinate
! ... (only for cylindrical coordinates)
!
            IF((itc == 0) .OR. (job_type == '3D')) THEN
              x(:) = 1.D0
              xb(:) = 1.D0
              inx(:) = 1.D0
              inxb(:) = 1.D0
            ELSE IF ((itc == 1) .AND. (job_type == '2D')) THEN
              xb(1)=0.D0
              x(1)=-0.5D0*dx(1)+xb(1)
              DO i=2,nx
                xb(i)=(xb(i-1)+dx(i))
                x(i)=(xb(i)-0.5D0*dx(i))
              END DO
              IF (xb(1) == 0.0D0) THEN
                inxb(1) = VERYBIG
              ELSE
                inxb(1) = 1.0D0/xb(1)
              END IF
              DO i=2,nx
                inxb(i)=1.D0/xb(i)
              END DO  
              DO i=1,nx
                inx(i)=1.D0/x(i)
              END DO  
            END IF
!
! ... Set "z" dimension
! ... (useful to set atmospheric conditions)
!
            zb(1) = zzero
            DO k=1,(nz-1)
              zb(k+1)=zb(k)+dz(k+1)
            END DO
!
! ... inverse of the cell sizes
!
            indx = 1.D0 / dx

            indz = 1.D0 / dz

            IF (job_type == '3D') THEN
              indy = 1.D0 / dy
            ELSE IF (job_type == '2D') THEN
              indy = 0.D0
            END IF
!
          RETURN 
        END SUBROUTINE grid_setup
!
!----------------------------------------------------------------------
!
      SUBROUTINE flic
!
        USE dimensions
        USE control_flags, ONLY: job_type
!
        IMPLICIT NONE
!
        INTEGER :: i, j, k, ij, n, ijk
!
        fl = 1
!
        IF( no <= 0 ) RETURN
!
! ...   no is the number of cell blocks (blunt body)
!
        IF( job_type == '2D' ) THEN

          DO n = 1, no
            DO k = iob(n)%zlo, iob(n)%zhi 
              DO i = iob(n)%xlo, iob(n)%xhi
                ijk = i + (k-1) * nx
                SELECT CASE ( iob(n)%typ )
                  CASE (2) 
                    fl(ijk) = 2
                  CASE (3) 
                    fl(ijk) = 3
                  CASE (4) 
                    fl(ijk) = 4
                  CASE (5) 
                    fl(ijk) = 5
                  CASE DEFAULT
                    fl(ijk) = 1
                END SELECT
              END DO
            END DO

          END DO

        ELSE IF( job_type == '3D' ) THEN

          DO n = 1, no
            DO k = iob(n)%zlo, iob(n)%zhi
              DO j = iob(n)%ylo, iob(n)%yhi
                DO i = iob(n)%xlo, iob(n)%xhi
                  ijk = i + ( (j-1) + (k-1)*ny ) * nx
                  SELECT CASE ( iob(n)%typ )
                    CASE (2)
                      fl(ijk) = 2
                    CASE (3)
                      fl(ijk) = 3
                    CASE (4)
                      fl(ijk) = 4
                    CASE (5)
                      fl(ijk) = 5
                    CASE DEFAULT
                      fl(ijk) = 1
                  END SELECT
                END DO
              END DO
            END DO
          END DO

        END IF

! ... MX3D CHECK THE CORNERS
!
        IF( job_type == '2D' ) THEN

          IF ( fl(nz*nx - 1) == 4 .AND. fl((nz-1)*nx) == 4) THEN
            fl(nz*nx) = 4
          END IF

        END IF
!
        RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
      END MODULE grid
!----------------------------------------------------------------------
