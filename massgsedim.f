!----------------------------------------------------------------------
      MODULE mass_ground
!----------------------------------------------------------------------
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
      USE domain_decomposition, ONLY: cell_owner, cell_g2l
      USE domain_mapping, ONLY: ncdom
      USE grid, ONLY: dx, dy, dz, r, xb, yb, zb, z
      USE kinds
      USE immersed_boundaries, ONLY: immb
      USE io_files, ONLY: tempunit, logunit
      USE io_parallel, ONLY: read_array
      USE postp_output, ONLY: improfile
      USE parallel, ONLY: mpime, root
      USE particles_constants, ONLY: rl
      USE postp_variables, ONLY: eps, time

      IMPLICIT NONE
!
      INTEGER :: i, j, k, ijk, imesh
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
      ALLOCATE(vf(ncdom))
      vf(:) = 1.D0
      IF (immb >= 1) THEN
        IF (mpime == root) &
          OPEN(tempunit,FILE='vf.dat',STATUS='OLD',FORM='UNFORMATTED')
        CALL read_array( tempunit, vf, dbl, .FALSE. )
        IF (mpime == root) CLOSE(tempunit)
      END IF
!
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
                imesh  = i + (j-1) * nx + (k-1) * nx * ny
                IF (cell_owner(imesh) == mpime) THEN
                  ijk = cell_g2l(imesh,mpime)
                  volume = dx(i) * dy(j) * dz(k) * vf(ijk)
                  mmap = mmap + SUM(eps(ijk,:)*rl(:)) * dz(k)
                  tmap(:) = tmap(:) + eps(ijk,:) * dz(k)
                  ! ... Map the value reached at any given position at
                  ! ... given time
              END IF
            END IF
          END DO
          m_map(i,j) = mmap
          t_map(i,j) = SUM(tmap(:))/0.67D0
          !
        END DO
      END DO
!
      CALL parallel_sum_real(m_map,nx*ny)
      CALL parallel_sum_real(t_map,nx*ny)
!
! ... Print out the map and the new 2D DEM file
!
      IF (mpime == root) THEN
        filnam1='m_map.'//lettera(tn)
        filnam2='t_map.'//lettera(tn)
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
      END IF
!
 122  FORMAT(10(1x,G14.6E3))

      DEALLOCATE(tmap)
      DEALLOCATE(m_map)
      DEALLOCATE(t_map)

      RETURN
      END SUBROUTINE massgs
!-----------------------------------------------------------------------
      END MODULE mass_ground
!----------------------------------------------------------------------
