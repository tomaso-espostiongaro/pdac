       SUBROUTINE STARTUP(NPROC,MPIME,ROOT,GROUP) 

       IMPLICIT NONE

#if defined MPI
       INCLUDE 'mpif.h'
#endif

       INTEGER NPROC
       INTEGER MPIME, I, ERR, ROOT, GROUP
!      ---------------------------
!      INITIALIZE MPI ENVIRONEMENT 
!      ---------------------------
       
#if defined MPI

       CALL MPI_INIT(ERR)  
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,ERR)
       CALL MPI_COMM_RANK(MPI_COMM_WORLD,MPIME,ERR)
       ROOT = 0
       GROUP = MPI_COMM_WORLD

#else

       NPROC = 1
       ROOT  = 0
       MPIME = 0
       GROUP = 0
        
#endif

       RETURN
       END
