!------------------------------------------------------------------------------!

      SUBROUTINE ISEND_REAL( array, n, idest, itag, ihand )
        IMPLICIT NONE
        REAL*8 array(*)
        INTEGER n, idest, itag, ihand
#if defined __MPI
        INCLUDE 'mpif.h'
        INTEGER :: ierr 
        CALL MPI_ISEND( array, n, MPI_DOUBLE_PRECISION, idest,
     &    itag, MPI_COMM_WORLD, ihand, ierr )
#endif
        RETURN
      END SUBROUTINE ISEND_REAL 


