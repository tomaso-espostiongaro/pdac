!----------------------------------------------------------------------
      MODULE sample_points
!----------------------------------------------------------------------
      USE postp_variables
!
      IMPLICIT NONE
      SAVE
!
! ... parameters for sampling
!
      INTEGER :: number_of_probes
      INTEGER, ALLOCATABLE :: ijk_probe(:), indx(:)
      LOGICAL :: assign_index
      CHARACTER(LEN=80) :: probe_file
!
      TYPE probe_point
            INTEGER :: nop
            INTEGER :: i
            INTEGER :: j
            INTEGER :: k
            REAL*8  :: x
            REAL*8  :: y
            REAL*8  :: z
       END TYPE probe_point
!
      TYPE(probe_point), ALLOCATABLE :: probe(:)
      INTEGER :: isamp
!
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE sample 
!
      USE control_flags, ONLY: job_type
      USE dimensions, ONLY: nsolid, ngas, nx, ny, nz
      USE domain_decomposition, ONLY: cell_owner, cell_g2l
      USE grid, ONLY: x, y, z
      USE io_files, ONLY: tempunit, logunit
      USE parallel, ONLY: mpime, root
!
      IMPLICIT NONE
      INTEGER :: ijk, imesh, i, j, k, n
      INTEGER :: nop, nv
      INTEGER, ALLOCATABLE :: imesh_probe(:), indx(:)
      CHARACTER(LEN = 11) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
      CHARACTER(LEN = 20) :: probenam
      !
      TYPE(probe_point), ALLOCATABLE :: probe(:)
!
! ... Allocate probes.
!
      ALLOCATE(probe(number_of_probes))
      ALLOCATE(imesh_probe(number_of_probes))
      ALLOCATE(indx(number_of_probes))
!
! ... Read probe file
! 
      IF (mpime == root) OPEN(tempunit, FILE=probe_file, STATUS='OLD',ERR=199)
      IF (mpime == root) READ(tempunit,*) 
      IF (mpime == root) WRITE(logunit,*) 'Sampling probes: '
      DO nop = 1, number_of_probes
        probe(nop)%nop = nop
        IF (assign_index) THEN
                IF (job_type == '2D') THEN
                        IF (mpime == root) READ(tempunit,*) probe(nop)%i,probe(nop)%k
                        CALL bcast_integer(probe(nop)%i,1,root)
                        CALL bcast_integer(probe(nop)%k,1,root)
                        i = probe(nop)%i
                        k = probe(nop)%k
                        imesh = i + (k-1)*nx
                        imesh_probe(nop) = imesh
                        !
                        probe(nop)%x = x(i)
                        probe(nop)%z = z(k)
                ELSE IF (job_type == '3D') THEN
                        IF (mpime == root) READ(tempunit,*) probe(nop)%i, probe(nop)%j, probe(nop)%k
                        CALL bcast_integer(probe(nop)%i,1,root)
                        CALL bcast_integer(probe(nop)%j,1,root)
                        CALL bcast_integer(probe(nop)%k,1,root)
                        i = probe(nop)%i
                        j = probe(nop)%j
                        k = probe(nop)%k
                        imesh = i + (j-1)*nx + (k-1)*nx*ny
                        imesh_probe(nop) = imesh
                        !
                        probe(nop)%x = x(i)
                        probe(nop)%y = y(j)
                        probe(nop)%z = z(k)
                END IF
        ELSE
                IF (job_type == '2D') THEN
                        IF (mpime == root) READ(tempunit,*) probe(nop)%x, probe(nop)%z
                        CALL bcast_real(probe(nop)%x,1,root)
                        CALL bcast_real(probe(nop)%z,1,root)
                        DO i=1,nx
                          IF (x(i) <= probe(nop)%x) probe(nop)%i = i
                        END DO
                        DO k=1,nz
                          IF (z(k) <= probe(nop)%z) probe(nop)%k = k
                        END DO
                        i = probe(nop)%i
                        k = probe(nop)%k
                        imesh = i + (k-1)*nx
                        imesh_probe(nop) = imesh
                ELSE IF (job_type == '3D') THEN
                        IF (mpime == root) READ(tempunit,*) probe(nop)%x, probe(nop)%y, probe(nop)%z
                        CALL bcast_real(probe(nop)%x,1,root)
                        CALL bcast_real(probe(nop)%y,1,root)
                        CALL bcast_real(probe(nop)%z,1,root)
                        DO i=1,nx
                          IF (x(i) <= probe(nop)%x) probe(nop)%i = i
                        END DO
                        DO j=1,ny
                          IF (y(j) <= probe(nop)%y) probe(nop)%j = j
                        END DO
                        DO k=1,nz
                          IF (z(k) <= probe(nop)%z) probe(nop)%k = k
                        END DO
                        i = probe(nop)%i
                        j = probe(nop)%j
                        k = probe(nop)%k
                        imesh = i + (j-1)*nx + (k-1)*nx*ny
                        imesh_probe(nop) = imesh
                END IF
        END IF
        IF (mpime == root .AND. job_type=='2D') THEN
          WRITE(logunit,'(3I4,4F8.2)') nop, i, k, x(i), z(k), probe(nop)%x, probe(nop)%z
        ELSE IF (mpime == root .AND. job_type=='3D') THEN
          WRITE(logunit,'(4I4,6F8.2)') nop, i, j, k, x(i), y(j), z(k), probe(nop)%x, probe(nop)%y, probe(nop)%z
        END IF
      END DO
      IF (mpime == root) CLOSE(tempunit)
!
! ... Sort the probes with progressively increasing index
! ... After the sorting, 'imesh_probe' contains the progressively
! ... increasing probe indexes, whereas 'indx' contains the
! ... probe indexes as read from the 'probe_file'
!
      CALL ikb07ad(imesh_probe(1), number_of_probes, indx(1))
!
      DO nop = 1, number_of_probes
        n = indx(nop)
        i = probe(n)%i
        j = probe(n)%j
        k = probe(n)%k
        IF (job_type == '3D') THEN
          probenam ='S'//lettera(n)//'_'//lettera(i)//'_'//&
                         lettera(j)//'_'//lettera(k)
        ELSE IF (job_type == '2D') THEN
          probenam ='S'//lettera(n)//'_'//lettera(i)//'_'//lettera(k)
        END IF
        OPEN(UNIT=tempunit, FILE=probenam, POSITION='APPEND')
        IF (cell_owner(imesh_probe(nop)) == mpime) THEN
          ijk = cell_g2l(imesh_probe(nop),mpime)
          IF (job_type == '3D') THEN
             WRITE(tempunit,100) time, p(ijk), &
             ug(ijk), vg(ijk), wg(ijk), tg(ijk), &
             (xgc(ijk,nv),nv=1, ngas), &
             (eps(ijk,nv), us(ijk,nv), & 
              vs(ijk,nv), ws(ijk,nv), &
              ts(ijk,nv), nv=1, nsolid), &
              rhom(ijk), mvm(ijk), pd(ijk)
          ELSE IF (job_type == '2D') THEN
             WRITE(tempunit,100) time, p(ijk), &
             ug(ijk), wg(ijk), tg(ijk), &
             (xgc(ijk,nv),nv=1, ngas), &
             (eps(ijk,nv), us(ijk,nv), &
              ws(ijk,nv), ts(ijk,nv),nv=1, nsolid), &
              rhom(ijk), mvm(ijk), pd(ijk)
          END IF
        END IF
        CLOSE(tempunit)
        CALL barrier
      END DO

 100  FORMAT( F8.2, 100(G14.6E3,1X) )
!
      DEALLOCATE(probe)
      DEALLOCATE(imesh_probe)
      DEALLOCATE(indx)
!
      RETURN
!
 199  CALL error ('sample','error in reading tempunit',tempunit)
!
      END SUBROUTINE sample
!----------------------------------------------------------------------
      END MODULE sample_points
!----------------------------------------------------------------------
