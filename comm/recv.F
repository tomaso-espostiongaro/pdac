!------------------------------------------------------------------------------!

      SUBROUTINE RECV_REAL( array, n, isour, itag )
        IMPLICIT NONE
        REAL*8 array(*)
        INTEGER n, isour, itag
#if defined __MPI
        INCLUDE 'mpif.h'
        INTEGER :: ierr 
        INTEGER :: istatus( MPI_STATUS_SIZE )
        CALL MPI_RECV( array, n, MPI_DOUBLE_PRECISION, isour,
     &    itag, MPI_COMM_WORLD, istatus, ierr )
#endif
        RETURN
      END SUBROUTINE RECV_REAL 


