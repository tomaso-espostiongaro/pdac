!------------------------------------------------------------------------------!
      SUBROUTINE MP_WALLTIME( walltime, mpime )
        IMPLICIT NONE
        REAL*8 :: walltime
        INTEGER :: mpime
#if defined __MPI
        INCLUDE 'mpif.h'
        INTEGER :: ierr 
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,MPIME,IERR)
        walltime = MPI_WTIME()
#else
        walltime = 0.0d0
        mpime = 0
#endif
        RETURN
      END SUBROUTINE MP_WALLTIME
!------------------------------------------------------------------------------!
