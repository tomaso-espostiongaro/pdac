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
      USE grid, ONLY: x, y, z
      USE io_files, ONLY: tempunit
!
      IMPLICIT NONE
      INTEGER :: ijk, i, j, k, n
      INTEGER :: nop, nv
      INTEGER, ALLOCATABLE :: ijk_probe(:), indx(:)
      CHARACTER(LEN = 11) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
      CHARACTER(LEN = 20) :: probenam
      !
      TYPE(probe_point), ALLOCATABLE :: probe(:)
!
! ... Allocate probes.
!
      ALLOCATE(probe(number_of_probes))
      ALLOCATE(ijk_probe(number_of_probes))
      ALLOCATE(indx(number_of_probes))
!
! ... Read probe file
! 
      OPEN(tempunit, FILE=probe_file, STATUS='OLD')
      DO nop = 1, number_of_probes
        probe(nop)%nop = nop
        IF (assign_index) THEN
                IF (job_type == '3D') THEN
                        READ(tempunit,*) probe(nop)%i, probe(nop)%j, probe(nop)%k
                        i = probe(nop)%i
                        j = probe(nop)%j
                        k = probe(nop)%k
                        ijk = i + (j-1)*nx + (k-1)*nx*ny
                        ijk_probe(nop) = ijk
                        !
                        probe(nop)%x = x(i)
                        probe(nop)%y = y(j)
                        probe(nop)%z = z(k)
                ELSE IF (job_type == '2D') THEN
                        READ(tempunit,*) probe(nop)%i,probe(nop)%k
                        i = probe(nop)%i
                        k = probe(nop)%k
                        ijk = i + (k-1)*nx
                        ijk_probe(nop) = ijk
                        !
                        probe(nop)%x = x(i)
                        probe(nop)%z = z(k)
                END IF
        ELSE
                IF (job_type == '3D') THEN
                        READ(tempunit,*) probe(nop)%x, probe(nop)%y, probe(nop)%z
                        probe(nop)%i = 1
                        probe(nop)%j = 1
                        probe(nop)%k = 1
                ELSE IF (job_type == '2D') THEN
                        READ(tempunit,*) probe(nop)%x, probe(nop)%z
                        probe(nop)%i = 1
                        probe(nop)%k = 1
                END IF
        END IF
      END DO
      CLOSE(tempunit)
!
! ... Sort the probes with progressively increasing index
! ... After the sorting, 'ijk_probe' contains the progressively
! ... increasing probe indexes, whereas 'indx' contains the
! ... probe indexes as read from the 'probe_file'
!
      CALL ikb07ad(ijk_probe(1), number_of_probes, indx(1))
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
        IF (job_type == '3D') THEN
           WRITE(tempunit,100) time, p(ijk_probe(nop)), &
           ug(ijk_probe(nop)), vg(ijk_probe(nop)), wg(ijk_probe(nop)), &
           tg(ijk_probe(nop)), &
           (xgc(ijk_probe(nop),nv),nv=1, ngas), &
           (eps(ijk_probe(nop),nv), us(ijk_probe(nop),nv), & 
            vs(ijk_probe(nop),nv), ws(ijk_probe(nop),nv), &
            ts(ijk_probe(nop),nv), nv=1, nsolid), &
            rm(ijk_probe(nop)), mvmp(ijk_probe(nop)), &
            pdp(ijk_probe(nop))
        ELSE IF (job_type == '2D') THEN
           WRITE(tempunit,100) time, p(ijk_probe(nop)), &
           ug(ijk_probe(nop)), wg(ijk_probe(nop)), &
           tg(ijk_probe(nop)), &
           (xgc(ijk_probe(nop),nv),nv=1, ngas), &
           (eps(ijk_probe(nop),nv), us(ijk_probe(nop),nv), &
            ws(ijk_probe(nop),nv), ts(ijk_probe(nop),nv),nv=1, nsolid), &
            rm(ijk_probe(nop)), mvmp(ijk_probe(nop)), &
            pdp(ijk_probe(nop))
        END IF
        CLOSE(tempunit)
      END DO

 100  FORMAT( F8.2, 100(G14.6E3,1X) )
!
      DEALLOCATE(probe)
      DEALLOCATE(ijk_probe)
      DEALLOCATE(indx)
!
      RETURN
      END SUBROUTINE sample
!----------------------------------------------------------------------
      END MODULE sample_points
!----------------------------------------------------------------------
