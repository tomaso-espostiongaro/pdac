!----------------------------------------------------------------------
      MODULE runtime_sampling
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
! ... parameters for sampling
!
      INTEGER :: number_of_probes
      INTEGER, ALLOCATABLE :: imesh_probe(:), ijk_probe(:)
      REAL*8, ALLOCATABLE :: pw(:,:)
!
      TYPE probe_point
            INTEGER :: nop
            INTEGER :: i
            INTEGER :: j
            INTEGER :: k
            REAL*8  :: x
            REAL*8  :: y
            REAL*8  :: z
            INTEGER :: owner
            INTEGER :: localindex
       END TYPE probe_point
!
      TYPE(probe_point), ALLOCATABLE :: probe(:)
      INTEGER :: isrt
!
      PRIVATE :: probe_point, probe, number_of_probes, ijk_probe
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE set_sampling
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions
      USE domain_decomposition, ONLY: cell_owner, cell_g2l
      USE grid, ONLY: x, y, z
      USE io_files, ONLY: tempunit, logunit, tempunit2
      USE parallel, ONLY: mpime, root
!
      IMPLICIT NONE
      INTEGER :: nop, i, j, k, ijk, imesh
      INTEGER, ALLOCATABLE :: indx(:)
!
! ... Read probe file
!
      IF (mpime == root) THEN
        OPEN(tempunit, FILE='probes.dat', STATUS='OLD')
        OPEN(tempunit2, FILE='probes_bis.dat', STATUS='UNKNOWN')
        READ(tempunit,*,ERR=199) number_of_probes
        WRITE(logunit,*) 'sampling mixture at probes: '
      END IF
      CALL bcast_integer(number_of_probes,1,root)
!
! ... Allocate probes.
!
      ALLOCATE(probe(number_of_probes))
      ALLOCATE(imesh_probe(number_of_probes))
      ALLOCATE(indx(number_of_probes))
      ALLOCATE(pw(number_of_probes,nsolid+ngas+6))
!
      DO nop = 1, number_of_probes
        probe(nop)%nop = nop
        !
        ! ... Read probe locations (either in coordinates or indices)
        IF (isrt == 1) THEN
          IF (mpime == root) READ(tempunit,*) probe(nop)%i, probe(nop)%j, probe(nop)%k
          CALL bcast_integer(probe(nop)%i,1,root)
          CALL bcast_integer(probe(nop)%j,1,root)
          CALL bcast_integer(probe(nop)%k,1,root)
          i = probe(nop)%i
          j = probe(nop)%j
          k = probe(nop)%k
          probe(nop)%x = x(i)
          probe(nop)%y = y(j)
          probe(nop)%z = z(k)
          IF (mpime == root) WRITE(tempunit2,*) probe(nop)%x, probe(nop)%y, probe(nop)%z
        ELSE IF (isrt >=2) THEN
          IF (mpime == root) READ(tempunit,*) probe(nop)%x, probe(nop)%y, probe(nop)%z
          CALL bcast_real(probe(nop)%x,1,root)
          CALL bcast_real(probe(nop)%y,1,root)
          CALL bcast_real(probe(nop)%z,1,root)
          DO i=1,nx
            IF (x(i)<=probe(nop)%x) probe(nop)%i=i
          END DO
          DO j=1,ny
            IF (y(j)<=probe(nop)%y) probe(nop)%j=j
          END DO
          DO k=1,nz
            IF (z(k)<=probe(nop)%z) probe(nop)%k=k
          END DO
          IF (mpime == root) WRITE(tempunit2,*) probe(nop)%i, probe(nop)%j, probe(nop)%k
        END IF
        imesh = probe(nop)%i + (probe(nop)%j-1)*nx + (probe(nop)%k-1)*nx*ny
        imesh_probe(nop) = imesh
      END DO
      IF (mpime == root) THEN
        CLOSE(tempunit2)
        CLOSE(tempunit)
      END IF
!
! ... Sort the probes with progressively increasing index
! ... After the sorting, 'imesh_probe' contains the progressively
! ... increasing probe indexes, whereas 'indx' contains the
! ... probe indexes as read from the 'probe_file'
!
!      CALL ikb07ad(imesh_probe(1), number_of_probes, indx(1))
!
      DO nop = 1, number_of_probes
        probe(nop)%owner = cell_owner(imesh_probe(nop))
        probe(nop)%localindex = cell_g2l(imesh_probe(nop),probe(nop)%owner)
        IF (mpime == root) &
          WRITE(logunit,*) 'Probe #', nop, 'Owner', probe(nop)%owner, &
                           'Local index', probe(nop)%localindex
      END DO
!
      RETURN
!
 199  CALL error('set_sampling','Error in reading probes file unit', tempunit)
! 
      END SUBROUTINE set_sampling
!----------------------------------------------------------------------
      SUBROUTINE sample_fields
!
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nsolid, ngas, nx, ny, nz
      USE mixture_fields, ONLY: rhom, um, vm, wm, tm, yd
      USE grid, ONLY: x, y, z
      USE io_files, ONLY: tempunit, testunit
      USE parallel, ONLY: mpime, root
      USE pressure_epsilon, ONLY: p
      USE time_parameters, ONLY: time
!
      IMPLICIT NONE
      INTEGER :: ijk, nop
      CHARACTER*4 :: lettera
      CHARACTER(LEN = 14) :: filnam
!
      DO nop = 1, number_of_probes
        IF (probe(nop)%owner == mpime) THEN
          ijk = probe(nop)%localindex
          pw(nop,1) = p(ijk)
          pw(nop,2) = rhom(ijk)
          pw(nop,3) = um(ijk)
          IF (job_type == JOB_TYPE_3D) THEN
            pw(nop,4) = vm(ijk)
          ELSE
            pw(nop,4) = 0.D0
          END IF
          pw(nop,5) = wm(ijk)
          pw(nop,6) = tm(ijk)
          pw(nop,7:) = yd(ijk,:)
        ELSE
          pw(nop,:) = 0.D0
        END IF
      END DO
      CALL parallel_sum_real(pw,number_of_probes*(nsolid+ngas+6))
!
      IF (mpime==root) THEN
        DO nop = 1, number_of_probes
          filnam = 'sample'//lettera(nop)//'.dat'
          OPEN(UNIT=tempunit, FILE=filnam, POSITION='APPEND')
          WRITE(tempunit,100) time, pw(nop,:), SUM(pw(nop,7:))
          CLOSE(tempunit)
        END DO
      END IF

 100  FORMAT( F10.5, 30(G14.6E3,1X) )
!
      RETURN
      END SUBROUTINE sample_fields
!----------------------------------------------------------------------
      END MODULE runtime_sampling
!----------------------------------------------------------------------
