!------------------------------------------------------------------------------!

      SUBROUTINE MP_WAIT( ihand )
        IMPLICIT NONE
        INTEGER ihand
#if defined MPI
        INCLUDE 'mpif.h'
        INTEGER :: ierr 
        INTEGER :: istatus( MPI_STATUS_SIZE )
        CALL MPI_WAIT( ihand, istatus, ierr )
#endif
        RETURN
      END SUBROUTINE MP_WAIT 


