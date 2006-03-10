!----------------------------------------------------------------------
      MODULE mass_ground
!----------------------------------------------------------------------
      USE postp_variables, ONLY: eps, time
      IMPLICIT NONE
!
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: m_map, t_map
      REAL*8   :: thickness
      INTEGER :: iground
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE massgs(tn)
! 
      USE control_flags, ONLY: job_type
      USE dimensions
      USE grid, ONLY: dx, dy, dz, r, xb, yb, zb, z
      USE immersed_boundaries, ONLY: immb
      USE io_files, ONLY: tempunit, logunit
      USE particles_constants, ONLY: rl
      USE io_serial, ONLY: improfile

      IMPLICIT NONE
!
      INTEGER :: i, j, k, ijk
      INTEGER, INTENT(IN) :: tn
      REAL*8  :: mmap, quota
      CHARACTER(LEN = 14) :: filnam1, filnam2
      CHARACTER(LEN = 4 ) :: lettera
      REAL*8 :: volume
      REAL*8, ALLOCATABLE :: vf(:)
      REAL*8, ALLOCATABLE :: tmap(:)

      IF (job_type == '2D') RETURN

      ALLOCATE(tmap(nsolid))
      ALLOCATE(m_map(nx,ny))
      ALLOCATE(t_map(nx,ny))
      tmap = 0.D0
      m_map = 0.D0
      t_map = 0.D0
!
      ALLOCATE(vf(ntot))
      vf(:) = 1.D0
      IF (immb >= 1) THEN
            OPEN(tempunit,FILE='vf.dat',STATUS='OLD',FORM='UNFORMATTED')
            READ(tempunit) (vf(ijk), ijk=1,ntot)
            CLOSE(tempunit)
      END IF
!
        filnam1='m_map.'//lettera(tn)
        filnam2='t_map.'//lettera(tn)
        WRITE(logunit,fmt="(/,'* Starting post-processing ',I5,' * ')" ) tn

        ! ... Compute the total solid mass to the ground
        !
        DO i = 1, nx
          DO j = 1, ny
            !
            tmap(:) = 0.D0
            mmap = 0.D0
            DO k = 1, nz
              quota = improfile(i,j,k)
              IF (quota >= 0.D0 .AND. quota<=thickness) THEN 
                  ijk  = i + (j-1) * nx + (k-1) * nx * ny
                  volume = dx(i) * dy(j) * dz(k) * vf(ijk)
                  mmap = mmap + SUM(eps(ijk,:)*rl(:)) * dz(k)
                  tmap(:) = tmap(:) + eps(ijk,:) * dz(k)
                  ! ... Map the value reached at any given position at
                  ! ... given time
              END IF
            END DO
            m_map(i,j) = mmap
            t_map(i,j) = SUM(tmap(:))/0.67D0
            !
          END DO
        END DO
        !
!

! ... Print out the map and the new 2D DEM file
!
      OPEN(UNIT=tempunit,FILE=filnam1)
      DO j = 1, ny
          WRITE(tempunit,122) (m_map(i,j), i=1, nx)
      END DO
      CLOSE(tempunit)
      OPEN(UNIT=tempunit,FILE=filnam2)
      DO j = 1, ny
          WRITE(tempunit,122) (t_map(i,j), i=1, nx)
      END DO
      CLOSE(tempunit)

 122  FORMAT(10(1x,G14.6E3))

      RETURN
      END SUBROUTINE massgs
!-----------------------------------------------------------------------
      END MODULE mass_ground
!----------------------------------------------------------------------
