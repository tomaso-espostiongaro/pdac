
      SUBROUTINE GATHER_COMPLEX(SEND,RECV,SIZE,root)
      IMPLICIT NONE
      INTEGER SIZE, root
      COMPLEX*16 SEND(SIZE)
      COMPLEX*16 RECV(SIZE)
#if defined MPI
      INCLUDE 'mpif.h'
      INTEGER IERR
      CALL MPI_GATHER(send,size,MPI_DOUBLE_COMPLEX,
     &   recv,size,MPI_DOUBLE_COMPLEX,root,MPI_COMM_WORLD,IERR)
#else
      RECV = SEND
#endif
      RETURN
      END

      SUBROUTINE GATHER_REAL(SEND,RECV,SIZE,root)
      IMPLICIT NONE
      INTEGER SIZE, root
      REAL*8 SEND(SIZE)
      REAL*8 RECV(SIZE)
#if defined MPI
      INCLUDE 'mpif.h'
      INTEGER IERR
      CALL MPI_GATHER(send,size,MPI_DOUBLE_PRECISION,
     &   recv,size,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,IERR)
#else
      RECV = SEND
#endif
      RETURN
      END

      SUBROUTINE GATHER_INTEGER(SEND,RECV,SIZE,root)
      IMPLICIT NONE
      INTEGER SIZE, root
      INTEGER SEND(SIZE)
      INTEGER RECV(SIZE)
#if defined MPI
      INCLUDE 'mpif.h'
      INTEGER IERR
      CALL MPI_GATHER(send,size,MPI_INTEGER,
     &   recv,size,MPI_INTEGER,root,MPI_COMM_WORLD,IERR)
#else
      RECV = SEND
#endif
      RETURN
      END
