!------------------------------------------------------------------------------!

      SUBROUTINE IRECV_REAL( array, n, isour, itag, ihand )
        IMPLICIT NONE
        REAL*8 array(*)
        INTEGER n, isour, itag, ihand
#if defined __MPI
        INCLUDE 'mpif.h'
        INTEGER :: ierr 
        CALL MPI_IRECV( array, n, MPI_DOUBLE_PRECISION, isour,
     &    itag, MPI_COMM_WORLD, ihand, ierr )
#endif
        RETURN
      END SUBROUTINE IRECV_REAL 


