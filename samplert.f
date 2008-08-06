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
      REAL*8, ALLOCATABLE :: pw(:)
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
      USE io_files, ONLY: tempunit, logunit
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
        READ(tempunit,*,ERR=199) number_of_probes
        WRITE(logunit,*) 'sampling pressure at probes: '
      END IF
      CALL bcast_integer(number_of_probes,1,root)
!
! ... Allocate probes.
!
      ALLOCATE(probe(number_of_probes))
      ALLOCATE(imesh_probe(number_of_probes))
      ALLOCATE(indx(number_of_probes))
      ALLOCATE(pw(number_of_probes))
!
      DO nop = 1, number_of_probes
        probe(nop)%nop = nop
        IF (job_type == JOB_TYPE_3D) THEN
                IF (mpime == root) THEN
                  READ(tempunit,*) probe(nop)%i, probe(nop)%j, probe(nop)%k
                  WRITE(logunit,*) probe(nop)%i, probe(nop)%j, probe(nop)%k, &
                      x(probe(nop)%i), y(probe(nop)%j), z(probe(nop)%k)
                END IF
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
        ELSE IF (job_type == JOB_TYPE_2D) THEN
                IF (mpime == root) THEN
                  READ(tempunit,*) probe(nop)%i,probe(nop)%k
                  WRITE(logunit,*) probe(nop)%i, probe(nop)%k
                END IF
                CALL bcast_integer(probe(nop)%i,1,root)
                CALL bcast_integer(probe(nop)%k,1,root)
                i = probe(nop)%i
                k = probe(nop)%k
                imesh = i + (k-1)*nx
                imesh_probe(nop) = imesh
                !
                probe(nop)%x = x(i)
                probe(nop)%z = z(k)
        END IF
      END DO
      IF (mpime == root) CLOSE(tempunit)
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
      END DO
!
      RETURN
!
 199  CALL error('set_sampling','error in reading temp unit', tempunit)
! 
      END SUBROUTINE set_sampling
!----------------------------------------------------------------------
      SUBROUTINE sample_pressure
!
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nsolid, ngas, nx, ny, nz
      USE grid, ONLY: x, y, z
      USE io_files, ONLY: tempunit, testunit
      USE parallel, ONLY: mpime, root
      USE pressure_epsilon, ONLY: p
      USE time_parameters, ONLY: time
!
      IMPLICIT NONE
      INTEGER :: ijk, nop
!
      DO nop = 1, number_of_probes
        IF (probe(nop)%owner == mpime) THEN
          ijk = probe(nop)%localindex
          pw(nop) = p(ijk)
        ELSE
          pw(nop) = 0.D0
        END IF
      END DO
      CALL parallel_sum_real(pw,number_of_probes)
!
      IF (mpime==root) THEN
        OPEN(UNIT=tempunit, FILE='pwav.dat', POSITION='APPEND')
        WRITE(tempunit,100) time, pw
        CLOSE(tempunit)
      END IF

 100  FORMAT( F10.5, 20(G14.6E3,1X) )
!
      RETURN
      END SUBROUTINE sample_pressure
!----------------------------------------------------------------------
      END MODULE runtime_sampling
!----------------------------------------------------------------------
