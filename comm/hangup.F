       SUBROUTINE HANGUP
       IMPLICIT NONE
#if defined __MPI
       INCLUDE 'mpif.h'
       INTEGER IERR
       CALL MPI_FINALIZE(IERR)
#endif
       RETURN
       END
