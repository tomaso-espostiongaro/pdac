      SUBROUTINE scatter_real(SEND,RECV,SIZE,root)
      IMPLICIT NONE
      INTEGER SIZE, ROOT
      REAL*8 SEND(SIZE), RECV(SIZE)
#if defined __MPI
      include 'mpif.h'
      INTEGER IERR
      CALL MPI_SCATTER(send,size,MPI_DOUBLE_PRECISION,
     &  recv,size,MPI_DOUBLE_PRECISION,ROOT,MPI_COMM_WORLD,IERR)
#else
      RECV = SEND
#endif
      RETURN
      END

      SUBROUTINE scatter_integer(SEND,RECV,SIZE,root)
      IMPLICIT NONE
      INTEGER SIZE, ROOT
      INTEGER SEND(SIZE), RECV(SIZE)
#if defined __MPI
      include 'mpif.h'
      INTEGER IERR
      CALL MPI_SCATTER(send,size,MPI_INTEGER,
     &  recv,size,MPI_INTEGER,ROOT,MPI_COMM_WORLD,IERR)
#else
      RECV = SEND
#endif
      RETURN
      END
