      SUBROUTINE error(a,b,n)

!  this routine prints an error message and terminates the program
!  ----------------------------------------------
!  END manual

      IMPLICIT NONE

#if defined MPI
      include 'mpif.h'
#endif

! ... declare subroutine arguments
      CHARACTER(LEN=*) a, b
      INTEGER n

      INTEGER :: nproc, mpime, ierr

! ... declare function

!  end of declarations
!  ----------------------------------------------

#if defined MPI
      CALL mpi_comm_size(mpi_comm_world,nproc,ierr)
      CALL mpi_comm_rank(mpi_comm_world,mpime,ierr)
#else
      MPIME = 0
      NPROC = 1
#endif

! ... print the error message
      WRITE (8,100) mpime, a, b, n

! ... terminate the program
      CALL myflush( 6 )  ! flush output streams

      IF( mpime == 0 ) THEN
        OPEN(UNIT=15, FILE='CRASH', POSITION='append', STATUS='unknown')
        WRITE (15,100) mpime, a, b, n
        CLOSE(UNIT=15)
      END IF

#if defined MPI
      CALL mpi_finalize( ierr )
      IF ( ierr/=0 ) THEN
        CALL mpi_abort(mpi_comm_world, ierr)
      END IF
#endif


100   FORMAT (/,' *** from PE : ',I3,'  *** in routine ',A, 
     &         /,' *** error msg. : ',A,' *** code ',I5, 
     &         /,' *** aborting ***', /)
      STOP 'CRASH'
      END SUBROUTINE
