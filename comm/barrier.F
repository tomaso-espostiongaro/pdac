!------------------------------------------------------------------------------!
      SUBROUTINE BARRIER
        IMPLICIT NONE
#if defined __MPI
        INCLUDE 'mpif.h'
        INTEGER :: ierr 
        CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
#endif
        RETURN
      END SUBROUTINE BARRIER 
!------------------------------------------------------------------------------!
