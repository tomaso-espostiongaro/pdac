!----------------------------------------------------------------------
      MODULE mass_partition
!----------------------------------------------------------------------
      USE postp_variables, ONLY: eps, time
!
      IMPLICIT NONE
!
      TYPE subdomain 
        INTEGER :: i1
        INTEGER :: i2
        INTEGER :: j1
        INTEGER :: j2
        INTEGER :: k1
        INTEGER :: k2
      END TYPE subdomain 
!
      INTEGER :: number_of_boxes, imassn
      CHARACTER(LEN=80) :: boxes_file
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE massn 
! 
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions
      USE domain_decomposition, ONLY: cell_owner, cell_g2l
      USE domain_mapping, ONLY: ncdom
      USE grid, ONLY: dx, dy, dz, r, xb, yb, zb, itc
      USE kinds
      USE immersed_boundaries, ONLY: immb
      USE io_files, ONLY: tempunit, logunit
      USE io_parallel, ONLY: read_array
      USE parallel, ONLY: mpime, root
      USE particles_constants, ONLY: rl

      IMPLICIT NONE
!
      INTEGER :: i1, i2, j1, j2, k1, k2
      INTEGER :: is, n
      INTEGER :: ijk, i, j, k, imesh
      REAL*8 :: volume, pi, twopi, totalmass
      REAL*8, ALLOCATABLE :: vf(:)
      REAL*8, ALLOCATABLE :: smass(:,:)
      TYPE( subdomain ), ALLOCATABLE :: box(:)
      REAL*8 :: x1, x2, y1, y2, z1, z2
      CHARACTER(LEN = 14) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
!
      ALLOCATE(vf(ncdom))
      vf(:) = 1.D0
      IF (immb >= 1) THEN
        IF (mpime == root) &
          OPEN(tempunit,FILE='vf.dat',STATUS='OLD',FORM='UNFORMATTED')
        CALL read_array( tempunit, vf, sgl, .FALSE. )
        IF (mpime == root) CLOSE(tempunit)
      END IF
!
      pi = 4.D0 * ATAN(1.D0)
      twopi = 2.D0 * pi
!
      ! ... Allocate boxes.
      ! ... Read the file containing the indexes or coordinates of the boxes
      !
      ALLOCATE(box(number_of_boxes))
      ALLOCATE(smass(number_of_boxes,nsolid))
      box(:)%i1 = 1
      box(:)%i2 = 1
      box(:)%j1 = 1
      box(:)%j2 = 1
      box(:)%k1 = 1
      box(:)%k2 = 1
!
      IF (mpime == root) OPEN(tempunit, FILE=boxes_file, STATUS='OLD',ERR=199)
!
      DO n = 1, number_of_boxes
                IF (job_type == JOB_TYPE_2D) THEN
                        IF (mpime == root) READ(tempunit,*) x1, x2, z1, z2
                        CALL bcast_real(x1,1,root)
                        CALL bcast_real(x2,1,root)
                        CALL bcast_real(z1,1,root)
                        CALL bcast_real(z2,1,root)
                        DO i=1,nx
                          IF (xb(i) < x1) box(n)%i1 = i
                          IF (xb(i) < x2) box(n)%i2 = i
                        END DO
                        DO k=1,nz
                          IF (zb(k) < z1) box(n)%k1 = k
                          IF (zb(k) < z2) box(n)%k2 = k
                        END DO
                        box(n)%j1 = 1
                        box(n)%j2 = 1
                ELSE IF (job_type == JOB_TYPE_3D) THEN
                        IF (mpime == root) READ(tempunit,*) x1, x2, y1, y2, z1, z2
                        CALL bcast_real(x1,1,root)
                        CALL bcast_real(x2,1,root)
                        CALL bcast_real(y1,1,root)
                        CALL bcast_real(y2,1,root)
                        CALL bcast_real(z1,1,root)
                        CALL bcast_real(z2,1,root)
                        DO i=1,nx
                          IF (xb(i) < x1) box(n)%i1 = i
                          IF (xb(i) < x2) box(n)%i2 = i
                        END DO
                        DO j=1,ny
                          IF (yb(j) < y1) box(n)%j1 = j
                          IF (yb(j) < y2) box(n)%j2 = j
                        END DO
                        DO k=1,nz
                          IF (zb(k) < z1) box(n)%k1 = k
                          IF (zb(k) < z2) box(n)%k2 = k
                        END DO
                END IF
      END DO
!
      IF (mpime == root) THEN
        CLOSE(tempunit)
        WRITE(logunit,*) 'Box limits'
        WRITE(logunit,*) (box(n)%i1, box(n)%i2, box(n)%j1, box(n)%j2, box(n)%k1,&
                    box(n)%k2, n=1,number_of_boxes)
      END IF
!
      ! ... Compute the total solid mass in a box
      !
      smass(:,:) = 0.D0
      totalmass = 0.D0
      DO n = 1, number_of_boxes
        i1 = box(n)%i1
        i2 = box(n)%i2
        j1 = box(n)%j1
        j2 = box(n)%j2
        k1 = box(n)%k1
        k2 = box(n)%k2
        DO k = k1, k2
          DO j = j1, j2
            DO i = i1, i2
              volume = 0.D0
              IF (job_type == JOB_TYPE_2D) THEN
                imesh = i + (k-1)*nx
                IF (cell_owner(imesh) == mpime) THEN
                  ijk = cell_g2l(imesh,mpime)
                  IF (itc == 1) THEN
                    volume = r(i) * dx(i) * dz(k) * vf(ijk)
                    volume = twopi * volume
                  ELSEIF (itc == 0) THEN
                    volume = dx(i) * dz(k)
                  END IF
                END IF
              ELSE IF (job_type == JOB_TYPE_3D) THEN
                imesh = i + (j-1)*nx + (k-1)*nx*ny
                IF (cell_owner(imesh) == mpime) THEN
                  ijk = cell_g2l(imesh,mpime)
                  volume = dx(i) * dy(j) * dz(k) * vf(ijk)
                END IF
              END IF
              IF (i/=1 .AND. i/=nx .AND. k/=1 .AND. k/=nz) THEN
                IF ((j/=1 .AND. j/=ny) .OR. (job_type == JOB_TYPE_2D)) THEN
                        DO is = 1, nsolid
                          smass(n,is)  = smass(n,is) + eps(ijk,is) * rl(is) * volume
                        END DO
                END IF
              END IF
            END DO
          END DO
        END DO
      END DO
!
      CALL parallel_sum_real(smass,number_of_boxes*nsolid)
      DO n = 1, number_of_boxes
        DO is = 1, nsolid
          totalmass = totalmass + smass(n,is)
        END DO
      END DO
!
      IF (mpime == root) THEN
        OPEN(tempunit, FILE='masspart.dat', STATUS='UNKNOWN', &
                                            POSITION='APPEND')
        WRITE(tempunit,100) &
          time, ((smass(n,is),is=1,nsolid),n=1,number_of_boxes),totalmass
        CLOSE(tempunit)
      END IF
 100  FORMAT(100(G30.15E3))

      RETURN
!
 199  CALL error ('massn','error in reading tempunit',tempunit)
!

      END SUBROUTINE massn
!-----------------------------------------------------------------------
      END MODULE mass_partition
!----------------------------------------------------------------------
