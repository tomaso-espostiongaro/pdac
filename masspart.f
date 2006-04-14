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
      USE dimensions
      USE grid, ONLY: dx, dy, dz, r, xb, yb, zb
      USE immersed_boundaries, ONLY: immb
      USE io_files, ONLY: tempunit
      USE particles_constants, ONLY: rl

      IMPLICIT NONE
!
      INTEGER :: i1, i2, j1, j2, k1, k2
      INTEGER :: is, ijk, i, j, k, n, counter
      REAL*8 :: volume, pi, twopi, totalmass
      REAL*8, ALLOCATABLE :: vf(:)
      REAL*8, ALLOCATABLE :: smass(:,:)
      TYPE( subdomain ), ALLOCATABLE :: box(:)
      REAL*8 :: x1, x2, y1, y2, z1, z2
      CHARACTER(LEN = 14) :: filnam
      CHARACTER(LEN = 4 ) :: lettera

      ALLOCATE(vf(ntot))
      vf(:) = 1.D0
      IF (immb >= 1) THEN
             OPEN(tempunit,FILE='vf.dat',STATUS='OLD',FORM='UNFORMATTED')
             READ(tempunit) (vf(ijk), ijk=1,ntot)
             CLOSE(tempunit)
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
      OPEN(tempunit, FILE=boxes_file, STATUS='OLD')
!
      DO n = 1, number_of_boxes
                IF (job_type == '2D') THEN
                        READ(tempunit,*) x1, x2, z1, z2
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
                ELSE IF (job_type == '3D') THEN
                        READ(tempunit,*) x1, x2, y1, y2, z1, z2
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
      CLOSE(tempunit)
      WRITE(*,*) 'Box limits'
      WRITE(*,*) (box(n)%i1, box(n)%i2, box(n)%j1, box(n)%j2, box(n)%k1,&
                  box(n)%k2, n=1,number_of_boxes)
!
      OPEN(tempunit, FILE='masspart.dat', STATUS='UNKNOWN', &
                                          POSITION='APPEND')
!      WRITE(tempunit,*) '****************************'

        ! ... Compute the total solid mass in a box
        !
        smass(:,:) = 0.D0
        totalmass = 0.D0
        DO n = 1, number_of_boxes
          counter = 0
          i1 = box(n)%i1
          i2 = box(n)%i2
          j1 = box(n)%j1
          j2 = box(n)%j2
          k1 = box(n)%k1
          k2 = box(n)%k2
          DO k = k1, k2
            DO j = j1, j2
              DO i = i1, i2
                IF (job_type == '2D') THEN
                  ijk = i + (k-1)*nx
                  volume = r(i) * dx(i) * dz(k) * vf(ijk)
                  volume = twopi * volume
                ELSE IF (job_type == '3D') THEN
                  ijk = i + (j-1)*nx + (k-1)*nx*ny
                  volume = dx(i) * dy(j) * dz(k) * vf(ijk)
                END IF
                IF (i/=1 .AND. i/=nx .AND. k/=1 .AND. k/=nz) THEN
                  IF ((j/=1 .AND. j/=ny) .OR. (job_type == '2D')) THEN
                          counter = counter + 1
                          DO is = 1, nsolid
                            smass(n,is)  = smass(n,is) + eps(ijk,is) * rl(is) * volume
                          END DO
                  END IF
                END IF
              END DO
            END DO
          END DO
          !WRITE(*,*) 'Block ', n, ' Counts: ', counter 
        END DO
        totalmass = SUM(smass(1:number_of_boxes,:))
        WRITE(tempunit,100) &
          time, ((smass(n,is),is=1,nsolid),n=1,number_of_boxes),totalmass
      CLOSE(tempunit)
 100  FORMAT(100(G30.15E3))

      RETURN
      END SUBROUTINE massn
!-----------------------------------------------------------------------
      END MODULE mass_partition
!----------------------------------------------------------------------
